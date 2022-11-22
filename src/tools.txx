
template<typename T>
int assert_vector_equality(T (&u)[3], T (&v)[3], T epsilon)
{
    int pass = 1;
    for (int i=0; i<3; i++)
    {
        if (std::abs(u[i]-v[i])>epsilon)
        {
            pass = 0;
        }
    }

    return pass;
}

template<typename T>
void cross(T (&u)[3], T (&v)[3], T (&w)[3])
{
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}