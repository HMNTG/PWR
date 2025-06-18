# Debugging a mpi program with n > 1 / mpiexec on windows
- generate visual studio solution

To set the number of MPI processes (the equivalent of `-n 2` in `mpiexec -n 2 mpi_program.exe`) directly in Visual Studio, follow these steps:

1.  **Right-click your project** in Solution Explorer and select **Properties**.
    
2.  Navigate to **Configuration Properties > Debugging**.
    
3.  In the **Command** field, enter the path to `mpiexec.exe` (e.g., `C:\Program Files\Microsoft MPI\Bin\mpiexec.exe`).
    
4.  In the **Command Arguments** field, enter:
    
    text
    
    `-n 2 "$(TargetPath)"`
    
    Replace `2` with the number of processes you want to start.
    

This tells Visual Studio to launch your program using `mpiexec` with the specified number of processes when you start debugging or running the project.

**Summary:**

- Go to project properties → Debugging
    
- Set Command to `mpiexec.exe`
    
- Set Command Arguments to `-n X "$(TargetPath)"` (replace X with your desired process count)[](https://learn.microsoft.com/en-us/answers/questions/736208/msmpi-for-parallel-computing-on-a-local-server-wit)
    

* * *

Visual Studio (since 2012) no longer includes built-in support for cluster/MPI debugging, so you cannot directly launch and debug all MPI processes from within the IDE as you would with a standard single-process application. However, you can still run your MPI program with more than one process (i.e., with `MPI_COMM_SIZE(MPI_COMM_WORLD) > 1`) by configuring the project to launch via `mpiexec` as described previously. This will correctly start multiple processes.

## Debugging MPI Programs in Visual Studio

**To debug an MPI program running with `mpiexec` in Visual Studio:**

1.  **Configure Project to Launch with mpiexec:**
    
    - Set the project’s Debugging Command to `mpiexec.exe` and Command Arguments to `-n 2 "$(TargetPath)"` (or however many processes you want).
2.  **Start Debugging:**
    
    - Start debugging (F5). This will launch multiple instances of your program as separate processes.
3.  **Attach the Debugger:**
    
    - In Visual Studio, go to **Debug > Attach to Process...**.
        
    - Find the running instances of your program (there will be one per MPI process).
        
    - Attach the debugger to one or more of these processes.
        
    
    > Tip: You may want to add a line like `getchar();` or a prompt for input at the beginning of your `main()` function, so the processes pause and give you time to attach the debugger before they proceed[5](http://supercomputingblog.com/mpi/debugging-an-mpi-application-with-microsoft-visual-studio/).
    
4.  **Resume Execution:**
    
    - Once attached, you can set breakpoints, step through code, and debug as usual. After attaching, continue execution from your prompt.

**Limitations:**

- You must manually attach to each process you want to debug.
    
- You cannot set breakpoints before process launch (unless using tricks like the input pause).
    
- Some extensions (like "ReAttach") can help speed up the attach process.
    

## find pid

- search for mpi / executable name
- take smallest PID, should be rank=0 process