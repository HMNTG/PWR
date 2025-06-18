# PWR25_1_Einleitung_Motivation

## S42 MPI-Programm zur parallelen Berechnung eines Skalarproduktes:

- [ ] completed

```c
int j, m, p, local_m, tag=0;
char message[100];
float local_dot, dot;
float local_x[100], local_y[100];
MPI_Status status;

MPI_Init(&argc, &argv);
MPI_Comm_rank( MPI_COMM_WORLD, &my_rank);
MPI_Comm_size( MPI_COMM_WORLD, &p);

local_m = m/p;
local_dot = 0.0;

for (j=0; j<local_m; j++)
{
    local_dot = local_dot + local_x[j] * local_y[j];
}

MPI_Reduce(&local_dot, &dot, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Finalize();
```

# PWR25_2_Datenverteilung_und_Matrix-Operationen

## S103 Parallele Matrix-Vektor Multiplikation: Zeilenverteilung

- [ ] completed

```c
local_n = n/p;
for (i=0; i<local_n; i++)
    local_c[i] = 0;
for (i=0; i<local_n; i++)
    for (j=0; j<m; j++)
        local_c[i] = local_c[i] + local_A[i][j] * b[j];
MPI_Allgather(local_c,local_n,MPI_DOUBLE,
              global_c,local_n,MPI_DOUBLE,comm);
```

- weiter varianten
    - Parallele Implementierung auf gemeinsamem Speicher
        - S106
    - Programmskizze: MV-Multiplikation, Linearkombinationen
        - S111

# PWR25_3_Loesung_linearer_Gleichungssysteme_mittels_Gauss-Elimination

## S133 Programmskizze zur Gauß-Elimination (sequentiell)

- [x] completed
- [ ] TODO need to abstract dynamic input

```c
double *gauss_sequential (double **a, double *b) {
    double *x, sum, l[MAX_SIZE];
    int i,j,k,r;
    x = (double *) malloc(n * sizeof(double));
    for (k = 0; k < n-1; k++) {
        r = max_col(a,k);
        if (k != r) exchange_row(a,b,r,k);
        for (i=(k+1); i < n; i++) {
            l[i] = a[i][k]/a[k][k];
            for (j=k; j < n; j++)
                a[i][j] = a[i][j] - l[i] * a[k][j];
            b[i] = b[i] - l[i] * b[k];
        }
    }
    for (k = n-1; k >= 0; k--)
    {
        sum = 0.0;
        for (j=k+1; j < n; j++)
            sum = sum + a[k][j] * x[j];
        x[k] = 1/a[k][k] * (b[k] - sum);
    }
    return x;
}
```

## S140 Programm zur Gauß-Elimination (zeilenzyklisch)

- [ ] completed

```c
double *gauss_cyclic (double **a, double *b)
{
    double *x, l[MAX_SIZE], *buf;
    int i,j,k,r, tag=42;
    MPI_Status status;
    struct { double val; int node; } z,y;

    x = (double *) malloc(n * sizeof(double));
    buf = (double *) malloc((n+1) * sizeof(double));
    for (k=0; k<n-1; k++) {
        r = max_col_loc(a,k);
        z.node = me;
        if (r != -1) z.val = abs(a[r][k]); else z.val = 0.0;
        MPI_Allreduce(&z,&y,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        if (k % p == y.node) { /* Pivotzeile und Zeile k auf einem Prozessor */
            if (k % p == me)
            {
                if (a[k][k] != y.val) exchange_row(a,b,r,k);
                copy_row(a,b,k,buf);
            }
        }
        else
    /* Pivotzeile und Zeile k auf unterschiedlichen Prozessoren */
    if (k % p == me) {
        copy_row(a,b,k,buf);
        MPI_Send(buf+k,n-k+1,MPI_DOUBLE,y.node,tag,
                 MPI_COMM_WORLD);
    }
    else if (y.node == me) {
        MPI_Recv(buf+k,n-k+1,MPI_DOUBLE,MPI_ANY_SOURCE,
                 tag,MPI_COMM_WORLD,&status);
        copy_exchange_row(a,b,r,buf);
    }
    MPI_Bcast(buf+k,n-k+1,MPI_DOUBLE,y.node,MPI_COMM_WORLD);
    if (((k % p) != y.node) && (k % p == me)) copy_row(a,b,buf,k);
    i = k+1; while (i % p != me) i++;
    for ( ; i<n; i+=p) {
        l[i] = a[i][k] / buf[k];
        for (j=k+1; j<n; j++)
            a[i][j] = a[i][j] - l[i]*buf[j];
        b[i] = b[i] - l[i]*buf[n];
    }
    for (k=n-1; k>=0; k--) {
        if (k % p == me) {
            sum = 0.0;
            for (j=k+1; j < n; j++) sum = sum + a[k][j] * x[j];
            x[k] = 1/a[k][k] * (b[k] - sum); }
        MPI_Bcast(&x[k],1,MPI_DOUBLE,k%p,MPI_COMM_WORLD);
    }
    return x;
}
```

## S157 Programmfragment zur Gauß-Elimination mit gesamtzyklischer Verteilung.

- [ ] completed

```c
double * gauss_double_cyclic (double **a, double *b)
{
    double *x, *buf, *elim_buf;
    int i,j,k,r,q, ql, size, buf_size, elim_size, psz;
    struct { double val; int pvtline; } z,y;
    MPI_Status status;

    x = (double *) malloc(n * sizeof(double));
    buf = (double *) malloc((n+1) * sizeof(double));
    elim_buf = (double *) malloc((n+1) * sizeof(double));
    for (k=0; k<n-1; k++) {
    if (member(me, Co(k))) {
        r = max_col_loc(a,k);
        z.pvtline = r; z.val = fabs(a[r][k]);
    }
    MPI_Allreduce(&z,&y,1,MPI_DOUBLE_INT,MPI_MAXLOC,comm(Co(k)));
    MPI_Bcast(&y,1,MPI_DOUBLE_INT,grp_leader(Co(k)),MPI_COMM_WORLD);
    r = y.pvtline;
    if (Ro(k) == Ro(r)) {
        /*Pivotzeile und Zeile k in einer Zeilengruppe */
        if (member(me, Ro(k))) {
            if (r != k) exchange_row_loc(a,b,r,k);
            copy_row_loc(a,b,k,buf); }
    }
    else /*Pivotzeile und Zeile k in unterschiedl. Zeilengruppen */
    if (member(me, Ro(k))) {
        copy_row_loc(a,b,buf);
        q = compute_partner(Ro(r),me);
        psz = compute_size(n,k,Ro(k));
        MPI_Send(buf+k,psz,MPI_DOUBLE,q,tag,MPI_COMM_WORLD); }
    else if (member(me,Ro(r))) {
        /* ausführender Prozessor enthält Teile der Pivotzeile */
        q = compute_partner(Ro(k),me);
        psz = compute_size(n,k,Ro(r));
        MPI_Recv(buf+k,psz,MPI_DOUBLE,q,tag,MPI_COMM_WORLD,&status);
        exchange_row_loc(a,b,r,buf);
    }
    for (q=0; q<p; q++) {
        if (member(q,Ro(r)) && member(me,Co(q))) {
            ql = rank(q,Cop(q)); buf_size = compute_size(n,k,Ro(k));
            MPI_Bcast(buf+k,buf_size,MPI_DOUBLE,ql,comm(Cop(q)));}
    }
    if ((Ro(k) != Ro(r)) && (member(me,Ro(k))))
        copy_row_loc(a,b,buf,k);
    if (member(me,Co(k))) elim_buf = compute_elim_fact_loc(a,b,k,buf);
    for (q=0; q<p; q++) {
        if (member(q,Co(k)) && member(me,Rop(q))) {
            ql = rank(q,Rop(q)); elim_size = compute_size(n,k,Co(k)));
            MPI_Bcast(elim_buf,elim_size,MPI_DOUBLE,ql,comm(Rop(q)));}
    }
    compute_local_entries(a,b,k,elim_buf,buf); }
    
    backward_substitution(a,b,x);
    
    return x;
}
```

# PWR25_4_Direkte_Verfahren_fuer_Gleichungssysteme_mit_Bandstruktur

- [ ] TODO kein code
- [ ] Gesamtalgorithmus zum Rekursiven Verdoppeln
    - [ ] alg overview S193
    - [ ] S195, S196 detail
- [ ] Zyklische Reduktion
    - [ ] Paralleler Algorithmus zur zyklischen Reduktion
        - [ ] alg S209 - S212
- [ ] Zyklische Reduktion für Bandmatrizen
    - [ ] Algorithmus der zyklischen Reduktion für Bandmatrizen
        - [ ] alg S223

# PWR25_5_Klassische_Iterationsverfahren_und_Cholesky-Faktorisierung

## S242 Jacobi-Verfahren

- [x] completed

```c
int Parallel_Jacobi(int n, int p, int max_it, float tol) {
    int i_loc, i_global, j, i, n_loc;
    double *x_tmp1[GLOB_MAX], *x_tmp2[GLOB_MAX];
    double *x_old, *x_new, *tmp;

    n_loc = n/p; /* lokale Blockgroesse */
    MPI_Allgather(loc_b, n_loc, MPI_DOUBLE, x_tmp1, n_loc,
                  MPI_DOUBLE, MPI_COMM_WORLD);

    x_new = x_tmp1; x_old = x_tmp2;

    do {
        tmp = x_new; x_new = x_old; x_old = tmp;
        for (i_loc = 0; i_loc < n_loc; i_loc++) {
            i_global = i_loc + me * n_loc;
            loc_x[i_loc] = loc_b[i_loc];

            for (j = 0; j < i_global; j++)
                loc_x[i_loc] = loc_x[i_loc] - loc_A[i_loc][j] * x_old[j];

            for (j = i_global+1; j < n; j++)
                loc_x[i_loc] = loc_x[i_loc] - loc_A[i_loc][j] * x_old[j];

            loc_x[i_loc] = loc_x[i_loc]/ loc_A[i_loc][i_global];
        }
        MPI_Allgather(loc_x, n_loc, MPI_DOUBLE, x_new, n_loc,
                      MPI_DOUBLE, MPI_COMM_WORLD);
    } while ((it_num < max_it) && (distance(x_old,x_new,n) >= tol));
}
```

## S244 Gauß-Seidel-Verfahren

- [x] completed

```c
n_local = n/p;
do {
    delta_x = 0.0;
    for (i = 0; i < n; i++) {
        s_k = 0.0;
        for (j = 0; j < n_local; j++)
            if (j + me * n_local != i)
                s_k = s_k + local_A[i][j] * x[j];
        root = i/n_local;
        i_local = i % n_local;
        MPI_Reduce(&s_k, &x[i_local], 1, MPI_FLOAT, MPI_SUM, root,
                   MPI_COMM_WORLD);
        if (me == root) {
            x_new = (b[i_local] - x[i_local]) / local_A[i][i_local];
            delta_x = max(delta_x, abs(x[i_local] - x_new));
            x[i_local] = x_new;
        }
    }
    MPI_Allreduce(&delta_x, &global_delta, 1, MPI_FLOAT,
                  MPI_MAX, MPI_COMM_WORLD);
} while(global_delta > tol);
```

## S253 Gauß-Seidel dünnbesetzt: parallele Implementierung

- [ ] completed

```c
sqn = sqrt(n);
do {
    for (l = 1; l <= sqn; l++) {
        for (j = me; j < l; j+=p) {
            i = l + j * (sqn-1) - 1; /* starte Numerierung bei 0 */
            x[i] = 0;
            if (i-sqn >= 0) x[i] = x[i] - a[i][i-sqn] * x[i-sqn];
            if (i > 0) x[i] = x[i] - a[i][i-1] * x[i-1];
            if (i+1 < n) x[i] = x[i] - a[i][i+1] * x[i+1];
            if (i+sqn < n) x[i] = x[i] - a[i][i+sqn] * x[i+sqn];
            x[i] = (x[i] + b[i]) / a[i][i]; }
        collect_elements(x,1); }
    for (l = 2; l <= sqn; l++) {
        for (j = me -1 +1; (j <= sqn -1) && (j >= 0); j+=p) {
            i = l * sqn + j * (sqn-1) - 1; x[i] = 0;
            if (i-sqn >= 0) x[i] = x[i] - a[i][i-sqn] * x[i-sqn];
            if (i > 0) x[i] = x[i] - a[i][i-1] * x[i-1];
            if (i+1 < n) x[i] = x[i] - a[i][i+1] * x[i+1];
            if (i+sqn < n) x[i] = x[i] - a[i][i+sqn] * x[i+sqn];
            x[i] = (x[i] + b[i]) / a[i][i]; }
        collect_elements(x,1); }
} while(convergence_test() < tol);
```

## S260 Rot-Schwarz-Anordnung

- [ ] completed

```c
local_nr = nr/p; local_ns = ns/p;
do {
    mestartr = me * local_nr;
    for (i= mestartr; i < mestartr + local_nr; i++) {
        xr[i] = 0;
        for (j ∈ N(i))
            xr[i] = xr[i] + a[i][j] * xs[j] ;
        xr[i] = (xr[i]+b[i]) / a[i][i] ;
    }
    collect_elements(xr);
    mestrarts = me * local_ns + nr;
    for (i= mestrarts; i < mestrarts + local_ns; i++) {
        xs[i] = 0;
        for (j ∈ N(i))
            xs[i] = xs[i] + a[i+nr][j] * xr[j];
        xs[i] = (xs[i] + b[i+nr]) / a[i+nr][i+nr];
    }
    collect_elements(xs);
} while (convergence_test());
```
