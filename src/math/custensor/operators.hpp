
#pragma once


namespace cut // Namespace Custom Tensor
{
    struct AddOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        { 
            return a + b;
        }
    };


    struct SubOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a - b;
        }
    };


    struct MultOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a * b;
        }
    };


    struct DivOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a, T b )
        {
            return a / b;
        }
    };


    struct ExpOp
    {
    public:
        // Define class methods
        template<typename T>
        static T apply( T a )
        {
            return std::exp( a );
        }
    };
}