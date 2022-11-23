
template<typename T>
inline int assert_vector_equality(int N, T* u, T* v, T epsilon)
{
    int pass = 1;
    for (int i=0; i<N; i++)
    {
        if (std::abs(u[i]-v[i])>epsilon)
        {
            pass = 0;
        }
    }

    return pass;
}


template<>
inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon)
{
    int pass = 1;
    for (int i=0; i<N; i++)
    {
        if ((u[i]-v[i]) != epsilon)
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


template<typename T>
void vAdd(int N, T* u, T* v, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = u[i] + v[i];
    }
}


template<typename T>
void vsAdd(int N, T* u, T s, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = s + u[i];
    }
}


template<typename T>
void vSub(int N, T* u, T* v, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = u[i] - v[i];
    }
}


template<typename T>
void vsSub(int N, T* u, T s, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = u[i] - s;
    }
}


template<typename T>
void vsSub(int N, T s, T* u, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = s - u[i];
    }
}