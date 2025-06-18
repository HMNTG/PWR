# PWR25_1_Einleitung_Motivation

## S42 MPI-Programm zur parallelen Berechnung eines Skalarproduktes:

- [ ] completed

```cpp
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

```cpp
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

```cpp
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

```cpp
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

```cpp
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