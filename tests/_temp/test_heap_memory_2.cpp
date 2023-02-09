

#include <iostream>

class P3;
void set_l1(P3* l1);
void set_l2(P3* l2);
void set_l3(P3* l3);


class P3
{
public:
    double a;
    double b;

    P3()
    {
        std::cout << "P3 - Constructor..." << std::endl;
    }

    ~P3()
    {
        std::cout << "P3 - Destructor..." << std::endl;
    }

};


class DatabaseStorage
{
public:
    P3* L1;
    P3* L2;
    P3* L3;

    DatabaseStorage()
    {
        std::cout << "DatabaseStorage - Constructor..." << std::endl;
    }

    ~DatabaseStorage()
    {
        std::cout << "DatabaseStorage - Destructor..." << std::endl;
        // delete L1;
        // delete L2;
        // delete L3;
    }

};


void check_databases(DatabaseStorage &dbs)
{
    std::cout << "Check database" << std::endl;
    std::cout << "L1->a: " << dbs.L1->a << std::endl;
    std::cout << "L1->b: " << dbs.L1->b << std::endl;
    std::cout << "L2->a: " << dbs.L2->a << std::endl;
    std::cout << "L2->b: " << dbs.L2->b << std::endl;
    std::cout << "L3->a: " << dbs.L3->a << std::endl;
    std::cout << "L3->b: " << dbs.L3->b << std::endl;
}


void set_databases(DatabaseStorage &dbs)
{
    std::cout << "Creating L1..." << std::endl;
    P3 L1 = P3();
    dbs.L1 = &L1;
    set_l1(dbs.L1);
    std::cout << "Creating L2..." << std::endl;
    P3 L2 = P3();
    dbs.L2 = &L2;
    set_l2(dbs.L2);
    std::cout << "Creating L2..." << std::endl;
    P3 L3 = P3();
    dbs.L3 = &L3;
    set_l3(dbs.L3);
    std::cout << "End building... " << std::endl;
    check_databases(dbs);
}


void set_l1(P3 *l1)
{
    l1->a = 0.0;
    l1->b = 1.0;
}


void set_l2(P3 *l2)
{
    l2->a = 2.0;
    l2->b = 3.0;
}


void set_l3(P3 *l3)
{
    l3->a = 4.0;
    l3->b = 5.0;
}


int main(void)
{
    std::cout << "DatabaseStorage - Creation...  " << std::endl;
    DatabaseStorage dbs = DatabaseStorage();
    std::cout << "DatbaseStorage - Created!" << std::endl;
    set_databases(dbs);
    std::cout << "Database Filled!" << std::endl;

    std::cout << "L1->a: " << dbs.L1->a << std::endl;
    std::cout << "L1->b: " << dbs.L1->b << std::endl;
    std::cout << "L2->a: " << dbs.L2->a << std::endl;
    std::cout << "L2->b: " << dbs.L2->b << std::endl;
    std::cout << "L3->a: " << dbs.L3->a << std::endl;
    std::cout << "L3->b: " << dbs.L3->b << std::endl;

    return 0;
}