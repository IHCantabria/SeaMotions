
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cstdlib>

// Include local modules
#include "../../src/tools.hpp"


class P3
{
public:
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    double e = 0.0;
    double f = 0.0;
    int num_vec = 100;
    double vec[100];

    P3()
    {
        // std::cout << "P3 - Constructor" << std::endl;
    }

    ~P3()
    {
        // std::cout << "P3 - Destructor" << std::endl;
    }
};


class IntegralsDb
{
public:
    bool is_heap = false;
    P3* l1;
    P3* l2;

    IntegralsDb()
    {
        // std::cout << "IntegralsDB - Constructor" << std::endl;
    }

    ~IntegralsDb()
    {
        // std::cout << "IntegralsDB - destructor"  << std::endl;
        if (this->is_heap)
        {
            // std::cout << " -> Deleting Heap memory - Destructor" << std::endl;
            delete this->l1;
            delete this->l2;
        }
    }
};


double operate_heap(IntegralsDb &dbs)
{
    double a = dbs.l1->a*dbs.l1->b;
    double b = dbs.l2->d*dbs.l2->c;
    double c = 0.0;
    for (int i=0; i<dbs.l1->num_vec; i++)
    {
        c += dbs.l1->vec[i]*dbs.l2->e*a*b;
    }

    return c;
}


void heap_memory(int num_iters)
{
    // Create storage object
    IntegralsDb dbs = IntegralsDb();
    dbs.l1 = new P3();
    dbs.l2 = new P3();
    dbs.is_heap = true;

    // Fill obejects
    dbs.l1->a = rand();
    dbs.l1->b = rand();
    dbs.l1->c = rand();
    dbs.l1->d = rand();
    dbs.l1->e = rand();
    dbs.l1->f = rand();

    dbs.l2->a = rand();
    dbs.l2->b = rand();
    dbs.l2->c = rand();
    dbs.l2->d = rand();
    dbs.l2->e = rand();
    dbs.l2->f = rand();

    // Loop over iterations performing operations
    double t0, t1;
    double elapsed_time_mean = 0.0;
    for (int i=0; i<num_iters; i++)
    {
        t0 = get_cpu_time();
        operate_heap(dbs);
        t1 = get_cpu_time();
        elapsed_time_mean += (t1-t0);
    }
    std::cout << "Elapsed time Heap: " << elapsed_time_mean/num_iters << std::endl;

}


void stack_memory(int num_iters)
{
    // Create storage object
    IntegralsDb dbs = IntegralsDb();
    P3 L1 = P3();
    dbs.l1 = &L1;
    P3 L2 = P3();
    dbs.l2 = &L2;
    dbs.is_heap = false;

    // Fill obejects
    dbs.l1->a = rand();
    dbs.l1->b = rand();
    dbs.l1->c = rand();
    dbs.l1->d = rand();
    dbs.l1->e = rand();
    dbs.l1->f = rand();

    dbs.l2->a = rand();
    dbs.l2->b = rand();
    dbs.l2->c = rand();
    dbs.l2->d = rand();
    dbs.l2->e = rand();
    dbs.l2->f = rand();

    // Loop over iterations performing operations
    double t0, t1;
    double elapsed_time_mean = 0.0;
    for (int i=0; i<num_iters; i++)
    {
        t0 = get_cpu_time();
        operate_heap(dbs);
        t1 = get_cpu_time();
        elapsed_time_mean += (t1-t0);
    }
    std::cout << "Elapsed time Stack: " << elapsed_time_mean/num_iters << std::endl;

}


int main(void)
{
    int N = 1e6;
    heap_memory(N);
    stack_memory(N);

    return 0;
}
