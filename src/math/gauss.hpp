
#ifndef __gauss_hpp
#define __gauss_hpp

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../config.hpp"
#include "math_tools.hpp"


/**************************************/
/********* Declare functions***********/
/**************************************/
void get_gauss_legendre(
                            int         num_points, 
                            cusfloat*   roots, 
                            cusfloat*   weights
                        );


/**************************************/
/********* Declare classes ************/
/**************************************/
struct GaussPoints
{
    // Define local variables
    int         np          = 0;
    cusfloat*   roots       = nullptr;
    cusfloat*   weights     = nullptr;

    // Define constructors and destructor
    GaussPoints( int np_in )
    {
        // Storage input data
        this->np        = np_in;

        // Allocate space for roots and weights
        this->roots     = generate_empty_vector<cusfloat>( this->np );
        this->weights   = generate_empty_vector<cusfloat>( this->np );

        // Get the required number of gauss points
        get_gauss_legendre(
                                this->np,
                                this->roots,
                                this->weights
                            );
        
    }

    ~GaussPoints( void )
    {
        mkl_free( this->roots );
        mkl_free( this->weights );
    }
    
};


template<std::size_t Order>
struct GaussPointsVec
{
    static_assert( Order >= 1 && Order < 20 );
};


template<>
struct GaussPointsVec<1>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 1;

    MEMALINGR static constexpr cusfloat       roots[1]          = { 0.0 };

    MEMALINGR static constexpr cusfloat       weights[1]        = { 2.0 };

};


template<>
struct GaussPointsVec<2>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 2;

    MEMALINGR static constexpr cusfloat       roots[2]          =   {
                                                                        -0.5773502691896257,
                                                                        0.5773502691896257
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[2]        =   { 
                                                                        1.0,
                                                                        1.0
                                                                    };

};


template<>
struct GaussPointsVec<3>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 3;

    MEMALINGR static constexpr cusfloat       roots[3]          =   {
                                                                        -0.7745966692414834,
                                                                        +0.0000000000000000,
                                                                        +0.7745966692414834
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[3]        =   { 
                                                                        0.5555555555555557,
                                                                        0.8888888888888883,
                                                                        0.5555555555555557
                                                                    };

};


template<>
struct GaussPointsVec<4>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 4;

    MEMALINGR static constexpr cusfloat       roots[4]          =   {
                                                                        -0.8611363115940526,
                                                                        -0.3399810435848563,
                                                                        +0.3399810435848563,
                                                                        +0.8611363115940526
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[4]        =   { 
                                                                        0.3478548451374538,
                                                                        0.6521451548625462,
                                                                        0.6521451548625462,
                                                                        0.3478548451374538
                                                                    };

};

template<>
struct GaussPointsVec<5>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 5;

    MEMALINGR static constexpr cusfloat       roots[5]          =   {
                                                                        -0.9061798459386640,
                                                                        -0.5384693101056831,
                                                                        +0.0000000000000000,
                                                                        +0.5384693101056831,
                                                                        +0.9061798459386640
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[5]        =   { 
                                                                        0.23692688505618897,
                                                                        0.47862867049936650,
                                                                        0.56888888888888910,
                                                                        0.47862867049936650,
                                                                        0.23692688505618897
                                                                    };

};

template<>
struct GaussPointsVec<6>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 6;

    MEMALINGR static constexpr cusfloat       roots[6]          =   {
                                                                        -0.93246951420315200,
                                                                        -0.66120938646626450,
                                                                        -0.23861918608319693,
                                                                        +0.23861918608319693,
                                                                        +0.66120938646626450,
                                                                        +0.93246951420315200
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[6]        =   { 
                                                                        0.17132449237917022,
                                                                        0.36076157304813855,
                                                                        0.46791393457269130,
                                                                        0.46791393457269130,
                                                                        0.36076157304813855,
                                                                        0.17132449237917022
                                                                    };

};


template<>
struct GaussPointsVec<7>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 7;

    MEMALINGR static constexpr cusfloat       roots[7]          =   {
                                                                        -0.9491079123427584,
                                                                        -0.7415311855993945,
                                                                        -0.4058451513773972,
                                                                        +0.0000000000000000,
                                                                        +0.4058451513773972,
                                                                        +0.7415311855993945,
                                                                        +0.9491079123427584
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[7]        =   { 
                                                                        0.12948496616886998,
                                                                        0.27970539148927664,
                                                                        0.38183005050511876,
                                                                        0.41795918367346910,
                                                                        0.38183005050511876,
                                                                        0.27970539148927664,
                                                                        0.12948496616886998
                                                                    };

};


template<>
struct GaussPointsVec<8>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 8;

    MEMALINGR static constexpr cusfloat       roots[8]          =   {
                                                                        -0.96028985649753630,
                                                                        -0.79666647741362670,
                                                                        -0.52553240991632900,
                                                                        -0.18343464249564984,
                                                                        +0.18343464249564984,
                                                                        +0.52553240991632900,
                                                                        +0.79666647741362670,
                                                                        +0.96028985649753630
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[8]        =   { 
                                                                        0.10122853629037560,
                                                                        0.22238103445337470,
                                                                        0.31370664587788766,
                                                                        0.36268378337836210,
                                                                        0.36268378337836210,
                                                                        0.31370664587788766,
                                                                        0.22238103445337470,
                                                                        0.10122853629037560
                                                                    };

};


template<>
struct GaussPointsVec<9>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 9;

    MEMALINGR static constexpr cusfloat       roots[9]          =   {
                                                                        -0.9681602395076261,
                                                                        -0.8360311073266358,
                                                                        -0.6133714327005904,
                                                                        -0.3242534234038089,
                                                                        +0.0000000000000000,
                                                                        +0.3242534234038089,
                                                                        +0.6133714327005904,
                                                                        +0.8360311073266358,
                                                                        +0.9681602395076261
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[9]        =   { 
                                                                        0.08127438836157415,
                                                                        0.18064816069485760,
                                                                        0.26061069640293550,
                                                                        0.31234707704000286,
                                                                        0.33023935500125967,
                                                                        0.31234707704000286,
                                                                        0.26061069640293550,
                                                                        0.18064816069485760,
                                                                        0.08127438836157415
                                                                    };

};


template<>
struct GaussPointsVec<10>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 10;

    MEMALINGR static constexpr cusfloat       roots[10]         =  {
                                                                        -0.97390652851717170,
                                                                        -0.86506336668898440,
                                                                        -0.67940956829902440,
                                                                        -0.43339539412924720,
                                                                        -0.14887433898163116,
                                                                        +0.14887433898163116,
                                                                        +0.43339539412924720,
                                                                        +0.67940956829902440,
                                                                        +0.86506336668898440,
                                                                        +0.97390652851717170
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[10]       =  { 
                                                                        0.06667134430868715,
                                                                        0.14945134915058056,
                                                                        0.21908636251598224,
                                                                        0.26926671930999674,
                                                                        0.29552422471475330,
                                                                        0.29552422471475330,
                                                                        0.26926671930999674,
                                                                        0.21908636251598224,
                                                                        0.14945134915058056,
                                                                        0.06667134430868715
                                                                    };

};


template<>
struct GaussPointsVec<11>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 11;

    MEMALINGR static constexpr cusfloat       roots[11]         =  {
                                                                        -0.97822865814605700,
                                                                        -0.88706259976809530,
                                                                        -0.73015200557404940,
                                                                        -0.51909612920681190,
                                                                        -0.26954315595234496,
                                                                        +0.00000000000000000,
                                                                        +0.26954315595234496,
                                                                        +0.51909612920681190,
                                                                        +0.73015200557404940,
                                                                        +0.88706259976809530,
                                                                        +0.97822865814605700
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[11]       =  { 
                                                                        0.055668567116174315,
                                                                        0.125580369464903860,
                                                                        0.186290210927734480,
                                                                        0.233193764591990500,
                                                                        0.262804544510246650,
                                                                        0.272925086777900600,
                                                                        0.262804544510246650,
                                                                        0.233193764591990500,
                                                                        0.186290210927734480,
                                                                        0.125580369464903860,
                                                                        0.055668567116174315
                                                                    };

};


template<>
struct GaussPointsVec<12>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 12;

    MEMALINGR static constexpr cusfloat       roots[12]         =  {
                                                                        -0.98156063424671920,
                                                                        -0.90411725637047490,
                                                                        -0.76990267419430470,
                                                                        -0.58731795428661750,
                                                                        -0.36783149899818013,
                                                                        -0.12523340851146897,
                                                                        +0.12523340851146897,
                                                                        +0.36783149899818013,
                                                                        +0.58731795428661750,
                                                                        +0.76990267419430470,
                                                                        +0.90411725637047490,
                                                                        +0.98156063424671920
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[12]       =  { 
                                                                        0.04717533638651318,
                                                                        0.10693932599531779,
                                                                        0.16007832854334605,
                                                                        0.20316742672306570,
                                                                        0.23349253653835467,
                                                                        0.24914704581340258,
                                                                        0.24914704581340258,
                                                                        0.23349253653835467,
                                                                        0.20316742672306570,
                                                                        0.16007832854334605,
                                                                        0.10693932599531779,
                                                                        0.04717533638651318
                                                                    };

};


template<>
struct GaussPointsVec<13>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 13;

    MEMALINGR static constexpr cusfloat       roots[13]         =  {
                                                                        -0.98418305471858810,
                                                                        -0.91759839922297790,
                                                                        -0.80157809073330990,
                                                                        -0.64234933944034030,
                                                                        -0.44849275103644687,
                                                                        -0.23045831595513483,
                                                                        +0.00000000000000000,
                                                                        +0.23045831595513483,
                                                                        +0.44849275103644687,
                                                                        +0.64234933944034030,
                                                                        +0.80157809073330990,
                                                                        +0.91759839922297790,
                                                                        +0.98418305471858810
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[13]       =  { 
                                                                        0.04048400476531615,
                                                                        0.09212149983772760,
                                                                        0.13887351021978760,
                                                                        0.17814598076194554,
                                                                        0.20781604753688862,
                                                                        0.22628318026289750,
                                                                        0.23255155323087406,
                                                                        0.22628318026289750,
                                                                        0.20781604753688862,
                                                                        0.17814598076194554,
                                                                        0.13887351021978760,
                                                                        0.09212149983772760,
                                                                        0.04048400476531615
                                                                    };

};

template<>
struct GaussPointsVec<14>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 14;

    MEMALINGR static constexpr cusfloat       roots[14]         =  {
                                                                        -0.98628380869681220,
                                                                        -0.92843488366357350,
                                                                        -0.82720131506976500,
                                                                        -0.68729290481168550,
                                                                        -0.51524863635815410,
                                                                        -0.31911236892788974,
                                                                        -0.10805494870734363,
                                                                        +0.10805494870734363,
                                                                        +0.31911236892788974,
                                                                        +0.51524863635815410,
                                                                        +0.68729290481168550,
                                                                        +0.82720131506976500,
                                                                        +0.92843488366357350,
                                                                        +0.98628380869681220
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[14]       =  { 
                                                                        0.035119460331752915,
                                                                        0.080158087159759540,
                                                                        0.121518570687902950,
                                                                        0.157203167158193630,
                                                                        0.185538397477937820,
                                                                        0.205198463721295520,
                                                                        0.215263853463157660,
                                                                        0.215263853463157660,
                                                                        0.205198463721295520,
                                                                        0.185538397477937820,
                                                                        0.157203167158193630,
                                                                        0.121518570687902950,
                                                                        0.080158087159759540,
                                                                        0.035119460331752915
                                                                    };

};

template<>
struct GaussPointsVec<15>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 15;

    MEMALINGR static constexpr cusfloat       roots[15]         =  {
                                                                        -0.98799251802048540
                                                                        -0.93727339240070600
                                                                        -0.84820658341042720
                                                                        -0.72441773136017010
                                                                        -0.57097217260853880
                                                                        -0.39415134707756340
                                                                        -0.20119409399743454
                                                                        +0.00000000000000000
                                                                        +0.20119409399743454
                                                                        +0.39415134707756340
                                                                        +0.57097217260853880
                                                                        +0.72441773136017010
                                                                        +0.84820658341042720
                                                                        +0.93727339240070600
                                                                        +0.98799251802048540
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[15]       =  { 
                                                                        0.03075324199611815,
                                                                        0.07036604748810715,
                                                                        0.10715922046717173,
                                                                        0.13957067792615430,
                                                                        0.16626920581699410,
                                                                        0.18616100001556220,
                                                                        0.19843148532711160,
                                                                        0.20257824192556134,
                                                                        0.19843148532711160,
                                                                        0.18616100001556220,
                                                                        0.16626920581699410,
                                                                        0.13957067792615430,
                                                                        0.10715922046717173,
                                                                        0.07036604748810715,
                                                                        0.03075324199611815
                                                                    };

};

template<>
struct GaussPointsVec<16>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 16;

    MEMALINGR static constexpr cusfloat       roots[16]         =  {
                                                                        -0.98940093499164990,
                                                                        -0.94457502307323260,
                                                                        -0.86563120238783180,
                                                                        -0.75540440835500300,
                                                                        -0.61787624440264380,
                                                                        -0.45801677765722740,
                                                                        -0.28160355077925890,
                                                                        -0.09501250983763745,
                                                                        +0.09501250983763745,
                                                                        +0.28160355077925890,
                                                                        +0.45801677765722740,
                                                                        +0.61787624440264380,
                                                                        +0.75540440835500300,
                                                                        +0.86563120238783180,
                                                                        +0.94457502307323260,
                                                                        +0.98940093499164990
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[16]       =  { 
                                                                        0.027152459411756466,
                                                                        0.062253523938647620,
                                                                        0.095158511682492310,
                                                                        0.124628971255533630,
                                                                        0.149595988816576380,
                                                                        0.169156519395002120,
                                                                        0.182603415044923300,
                                                                        0.189450610455068090,
                                                                        0.189450610455068090,
                                                                        0.182603415044923300,
                                                                        0.169156519395002120,
                                                                        0.149595988816576380,
                                                                        0.124628971255533630,
                                                                        0.095158511682492310,
                                                                        0.062253523938647620,
                                                                        0.027152459411756466
                                                                    };

};

template<>
struct GaussPointsVec<17>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 17;

    MEMALINGR static constexpr cusfloat       roots[17]         =  {
                                                                        -0.99057547531441740,
                                                                        -0.95067552176876780,
                                                                        -0.88023915372698600,
                                                                        -0.78151400389680140,
                                                                        -0.65767115921669080,
                                                                        -0.51269053708647690,
                                                                        -0.35123176345387636,
                                                                        -0.17848418149584777,
                                                                        +0.00000000000000000,
                                                                        +0.17848418149584777,
                                                                        +0.35123176345387636,
                                                                        +0.51269053708647690,
                                                                        +0.65767115921669080,
                                                                        +0.78151400389680140,
                                                                        +0.88023915372698600,
                                                                        +0.95067552176876780,
                                                                        +0.99057547531441740
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[17]       =  { 
                                                                        0.02414830286854719,
                                                                        0.05545952937398706,
                                                                        0.08503614831717916,
                                                                        0.11188384719340384,
                                                                        0.13513636846852553,
                                                                        0.15404576107681064,
                                                                        0.16800410215645020,
                                                                        0.17656270536699292,
                                                                        0.17944647035620678,
                                                                        0.17656270536699292,
                                                                        0.16800410215645020,
                                                                        0.15404576107681064,
                                                                        0.13513636846852553,
                                                                        0.11188384719340384,
                                                                        0.08503614831717916,
                                                                        0.05545952937398706,
                                                                        0.02414830286854719
                                                                    };

};

template<>
struct GaussPointsVec<18>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 18;

    MEMALINGR static constexpr cusfloat       roots[18]         =  {
                                                                        -0.99156516842093090,
                                                                        -0.95582394957139760,
                                                                        -0.89260246649755580,
                                                                        -0.80370495897252310,
                                                                        -0.69168704306035320,
                                                                        -0.55977083107394750,
                                                                        -0.41175116146284270,
                                                                        -0.25188622569150560,
                                                                        -0.08477501304173532,
                                                                        +0.08477501304173532,
                                                                        +0.25188622569150560,
                                                                        +0.41175116146284270,
                                                                        +0.55977083107394750,
                                                                        +0.69168704306035320,
                                                                        +0.80370495897252310,
                                                                        +0.89260246649755580,
                                                                        +0.95582394957139760,
                                                                        +0.99156516842093090
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[18]       =  { 
                                                                        0.02161601352648239,
                                                                        0.04971454889497083,
                                                                        0.07642573025488873,
                                                                        0.10094204410628753,
                                                                        0.12255520671147814,
                                                                        0.14064291467065074,
                                                                        0.15468467512626527,
                                                                        0.16427648374583267,
                                                                        0.16914238296314368,
                                                                        0.16914238296314368,
                                                                        0.16427648374583267,
                                                                        0.15468467512626527,
                                                                        0.14064291467065074,
                                                                        0.12255520671147814,
                                                                        0.10094204410628753,
                                                                        0.07642573025488873,
                                                                        0.04971454889497083,
                                                                        0.02161601352648239
                                                                    };

};


template<>
struct GaussPointsVec<19>
{
public:
    // Define public attributes
    MEMALINGR static constexpr std::size_t    N                 = 19;

    MEMALINGR static constexpr cusfloat       roots[19]         =  {
                                                                        -0.99240684384358440,
                                                                        -0.96020815213483000,
                                                                        -0.90315590361481780,
                                                                        -0.82271465653714280,
                                                                        -0.72096617733522940,
                                                                        -0.60054530466168100,
                                                                        -0.46457074137596090,
                                                                        -0.31656409996362990,
                                                                        -0.16035864564022534,
                                                                        +0.00000000000000000,
                                                                        +0.16035864564022534,
                                                                        +0.31656409996362990,
                                                                        +0.46457074137596090,
                                                                        +0.60054530466168100,
                                                                        +0.72096617733522940,
                                                                        +0.82271465653714280,
                                                                        +0.90315590361481780,
                                                                        +0.96020815213483000,
                                                                        +0.99240684384358440
                                                                    };

    MEMALINGR static constexpr cusfloat       weights[19]       =  { 
                                                                        0.019461788229728976,
                                                                        0.044814226765699290,
                                                                        0.069044542737641180,
                                                                        0.091490021622449710,
                                                                        0.111566645547333750,
                                                                        0.128753962539335850,
                                                                        0.142606702173606380,
                                                                        0.152766042065859230,
                                                                        0.158968843393954030,
                                                                        0.161054449848783250,
                                                                        0.158968843393954030,
                                                                        0.152766042065859230,
                                                                        0.142606702173606380,
                                                                        0.128753962539335850,
                                                                        0.111566645547333750,
                                                                        0.091490021622449710,
                                                                        0.069044542737641180,
                                                                        0.044814226765699290,
                                                                        0.019461788229728976
                                                                    };

};

#endif