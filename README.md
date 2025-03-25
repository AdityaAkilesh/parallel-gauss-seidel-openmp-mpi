# Parallel Gauss-Seidel Method using OpenMP and MPI

## 📌 Overview
This project implements the Gauss-Seidel iterative method to solve linear systems using **parallelization techniques**:
- **OpenMP** for shared memory systems (multi-threading)
- **MPI** for distributed memory systems (multi-node)

It compares execution time and convergence behavior between serial, OpenMP, and MPI implementations. All experiments were conducted on the **Nova Cluster** using the **GCC compiler toolchain** for compiling both OpenMP and MPI code.

---

## 📁 Project Structure
```
parallel-gauss-seidel-openmp-mpi/
|
├── Gauss Siedel OMP/             # OpenMP Implementation
│   ├── gs_omp.c                  # OpenMP-based solver
│   ├── gs_omp.ipynb              # Analysis notebook
│   ├── data.csv                  # Timing results (serial vs parallel)
│   ├── outputX.txt               # Output files per matrix size/thread
│   ├── Makefile, myjob           # GCC Makefile & SLURM job submission script
│   └── *.png                     # Plots for visualization
│
├── Gauss Siedel MPI/             # MPI Implementation
│   ├── gs_mpi.c                  # MPI-based solver
│   ├── gs_mpi.ipynb              # Jupyter notebook (MPI result analysis)
│   ├── data.csv                  # Timing results (processes vs size)
│   ├── outputX.data              # Output per test run
│   ├── Makefile, myjob           # GCC Makefile & SLURM job submission script
│   └── *.png                     # Results and comparison plots
|
├── results.png                   # Combined summary result chart
├── LICENSE                       # Open-source license (optional)
├── CONTRIBUTING.md               # Contribution guidelines
└── README.md                     # Documentation
```

---

## ⚙️ Build & Execution (on Nova Cluster)
All programs were compiled using **GCC** and executed via **SLURM batch jobs** on the Nova HPC cluster.

### 🧵 OpenMP
Build using the provided Makefile:
```bash
make omp      # builds main.exe for OpenMP using gcc -fopenmp
```
Submit via SLURM:
```bash
sbatch myjob   # Adjust job script with thread and matrix size parameters
```
Manual run:
```bash
./main.exe <thread_count> <matrix_size>
```
Example:
```bash
./main.exe 4 512
```

### 💻 MPI
Build using:
```bash
make mpi      # builds main.exe using mpicc (GCC wrapper)
```
Submit via SLURM:
```bash
sbatch myjob   # Update process count and matrix size as needed
```
Or run interactively:
```bash
mpirun -np <num_processes> ./main.exe <matrix_size>
```
Example:
```bash
mpirun -np 4 ./main.exe 512
```

---

## 📝 Sample SLURM Job Script (myjob)
```bash
#!/bin/bash
#SBATCH --job-name=gs_parallel
#SBATCH --output=output.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --partition=short

module load gcc
module load openmpi

mpirun ./main.exe 512
```
Modify `--ntasks` and matrix size as needed. For OpenMP jobs, replace `mpirun` with `./main.exe <threads> <size>` and load OpenMP modules if needed.

---

## 📊 Data Format
The implementations append timing results to `data.csv`:

### OpenMP:
```
thread_count;MatrixSize;Total_time_serial;total_time_parallel
4;128;0.013517;0.004629
```

### MPI:
```
Processes;Size;Total-time;
4;512;0.007229;
```

---

## 📈 Visual Analysis
- Timing data is processed using Python notebooks and visualized with bar/line plots
- Comparisons are made across matrix sizes, thread/process counts, and convergence iterations

![](Gauss%20Siedel%20OMP/gs_omp%20result.png)
![](Gauss%20Siedel%20MPI/results.png)

---

## 🧠 Implementation Details

### ✅ OpenMP Highlights (Red-Black Ordering)
- Matrix updated in two passes: red (odd sum of indices) and black (even sum)
- `#pragma omp parallel for` with reduction on convergence difference
- Timing captured using `omp_get_wtime()`
- Output written to `data.csv`

### ✅ MPI Highlights
- Master initializes full matrix, distributes blocks using `MPI_Scatterv`
- Each process solves its block independently with Gauss-Seidel
- Results gathered using `MPI_Gatherv`
- Execution time measured using `MPI_Wtime()`

---




