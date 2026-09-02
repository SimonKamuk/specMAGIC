# OpenMP quickstart

This code uses **OpenMP** to parallelize selected computationally expensive loops. OpenMP is a compiler-supported API for shared-memory parallel programming in C and C++.

OpenMP is a lightweight way of distributing work across multiple CPU threads without manually creating and managing thread objects.

It may be useful to see the Quickstart section of the page on C++ programming. A more thorough guide to OpenMP, including examples and compilation, is [provided here](https://github-pages.ucl.ac.uk/research-computing-with-cpp/08openmp/02_intro_openmp.html), and [these](https://www.openmp.org/wp-content/uploads/omp-hands-on-SC08.pdf) are also useful slides.

## Compiler pragmas

OpenMP is primarily controlled through pragmas.

A pragma is a compiler directive that begins with `#pragma`. Pragmas provide instructions to the compiler but are not part of the C++ language itself.

```cpp
#pragma omp parallel 
{
    // code here is executed by multiple threads 
}
```

The OpenMP runtime recognises pragmas beginning with `omp` and generate the appropriate parallel code.

If OpenMP support is disabled during compilation, these pragmas are ignored and the code executes serially.

## Setting the number of threads

The number of threads used by OpenMP can be controlled with:

```cpp
omp_set_num_threads(n_threads);

// Example
omp_set_num_threads(8);
```

This requests that future OpenMP parallel regions use up to 8 threads.

Typically this is called once during initialisation.

## Querying the number of threads

Inside a parallel region, the number of active threads can be obtained using `omp_get_num_threads()` . The thread ID can be gotten with `omp_get_thread_num()`.

```cpp
omp_set_num_threads(4);
#pragma omp parallel
{
    int n_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
    printf("Hello from thread %d of %d! \n", thread_id, n_threads);
} 
>> Hello from thread 2 of 4!
>> Hello from thread 1 of 4!
>> Hello from thread 0 of 4!
>> Hello from thread 3 of 4!
```

The order in which the threads write from within a parallel region is not deterministic.

## Parallel Regions

A parallel region is created using

```cpp
#pragma omp parallel
{
    // executed by all threads
}
```

When execution reaches this block:

1. OpenMP creates a team of worker threads.
2. Each thread executes the code inside the block.
3. The threads synchronize at the end of the block.
4. Execution continues with a single thread.

## Parallelising loops

The most common OpenMP construct is the parallelized `for` loop:

```cpp
#pragma omp parallel for 
for (int i = 0; i < N; i++) {     
    doStuff(i); 
} 
```

The iterations are automatically divided among the available threads.

For example, with 4 threads and 100 iterations, each thread may process approximately 25 iterations.

## When `parallel` and `for` are separate

OpenMP provides a convenient combined form of `parallel for`, see above.

However, this project intentionally uses separate `parallel` and `for` directives:

```cpp
#pragma omp parallel
{
    // Thread-local setup

    #pragma omp for
    for (int i = 0; i < N; i++)
    {
        doStuff(i);
    }

    // Additional parallel work
}
```

This allows scratch buffers to be allocated once per thread and reused across multiple loops, reducing allocation overhead and improving performance.

For example:

```cpp
#pragma omp parallel
{
    Workspace workspace;

    #pragma omp for
    for (int i = 0; i < N; ++i)
    {
        compute_a(i, workspace);
    }

    #pragma omp for
    for (int i = 0; i < N; ++i)
    {
        compute_b(i, workspace);
    }
}
```

Each thread keeps its own `workspace` object for the entire lifetime of the parallel region. This is particularly advantageous if there are many more iterations than there are threads, as it reduces memory overhead significantly.

## Shared vs thread-local variables

OpenMP uses a **s**hared-memory model.

All threads can access the same process memory, which means care must be taken when writing data.

### Safe: independent writes

```cpp
#pragma omp for
for (int i = 0; i < N; i++)
{
    output[i] = compute(i);
}
```

Each iteration writes to a different element of `output`, so no threads interfere with each other.

### Unsafe: shared writes

```cpp
double total = 0.0;

#pragma omp for
for (int i = 0; i < N; ++i)
{
    total += values[i];
}
```

Multiple threads may attempt to update `total` simultaneously, producing incorrect results.

This is known as a race condition.

Any variable that is being written to should generally be **thread-local** unless explicit synchronization is used.

Good:

```cpp
#pragma omp parallel
{
    Workspace workspace;  // private to each thread

    ...
}
```

Potentially problematic:

```cpp
Workspace workspace;

#pragma omp parallel
{
    workspace.update(...);
}
```

In the second example, multiple threads modify the same object concurrently.

When designing parallel code:

* Read-only data, e.g. look-up tables, may safely be shared. This does not degrade performance.
* Data written by a thread should usually be owned exclusively by that thread.
* Shared writable data requires synchronisation, which should be avoided unless necessary.