#ifndef SIM_VECTOR_C
#define SIM_VECTOR_C 

#include "immintrin.h"
//TODO: these need to be global between files at least refactor N and size to be passed directly
#define N 200
#define SIZE (N + 2) * (N + 2)
#define IX(i,j) ((j) + (N+2)*(i))
#define SWAP(x0,x) {float *tmp = x0; x0=x;x=tmp;}
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

//Todo: vectorize this, just need to handle a mask for the extra values
void add_source(int n, float* x, float* s, float dt){
    int size = (n+2) * (n+2);
    for(int i = 0; i < size; i++){
        x[i] += dt * s[i];  
    }
}

void set_bnd(int n, int b, float* x){
    for(int i = 1; i < n; i++){
        x[IX(0  ,i)] = ((b==1) * -x[IX(1,i)]) + ((b!=1) * x[IX(1,i)]); 
        x[IX(n+1,i)] = ((b==1) * -x[IX(n,i)]) + ((b!=1) * x[IX(n,i)]); 
        x[IX(i,0  )] = ((b==2) * -x[IX(i,1)]) + ((b!=2) * x[IX(i,1)]);
        x[IX(i,n+1)] = ((b==2) * -x[IX(i,n)]) + ((b!=2) * x[IX(i,n)]); 
        //x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
        //x[IX(n+1,i)] = b==1 ? -x[IX(n,i)] : x[IX(n,i)];
        //x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        //x[IX(i,n+1)] = b==2 ? -x[IX(i,n)] : x[IX(i,n)];
    }
    x[IX(0  ,0  )] = 0.5 * (x[IX(1,0  )] + x[IX(0  ,1)]);
    x[IX(0  ,n+1)] = 0.5 * (x[IX(1,n+1)] + x[IX(0  ,n)]);
    x[IX(n+1,0  )] = 0.5 * (x[IX(n,0  )] + x[IX(n+1,1)]);
    x[IX(n+1,n+1)] = 0.5 * (x[IX(n,n+1)] + x[IX(n+1,n)]);
}



//Seems it's no better than 256, sounds like the instructions are still encoded as 2 avx2 calls
void diffuse_512(int n, int b, float* x, float* x0, float diff, float dt){
    float a = dt * diff * n * n;
    float a_div = 1/(1+4 * a);
    float tmp[SIZE];
    __m512 _top, _bottom, _left, _right, _resA, _resB, _a, _a_div, _center;
    _a = _mm512_set1_ps(a);
    _a_div = _mm512_set1_ps(a_div);
    for(int k = 0; k < 20; k++){
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j+=16){
                _top = _mm512_loadu_ps(&x[IX(i-1, j)]);
                _bottom = _mm512_loadu_ps(&x[IX(i+1, j)]);
                _left = _mm512_loadu_ps(&x[IX(i, j-1)]);
                _right = _mm512_loadu_ps(&x[IX(i, j+1)]);
                _center = _mm512_loadu_ps(&x0[IX(i, j)]);
                _resA = _mm512_add_ps(_top,_bottom);
                _resB = _mm512_add_ps(_left,_right);
                _resA = _mm512_add_ps(_resA,_resB);
                _resA = _mm512_fmadd_ps(_resA, _a , _center);
                _resA = _mm512_mul_ps(_resA, _a_div);

                _mm512_storeu_ps(&tmp[IX(i,j)], _resA);
            }
        }
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j++){
                x[IX(i,j)] = tmp[IX(i,j)];
            }
        }
        set_bnd(n, b, x);
    }
}

//can I avx the i at all, or better memoize
//#pragma optimize("", off)
void diffuse(int n, int b, float* x, float* x0, float diff, float dt){
    float a = dt * diff * n * n;
    float a_div = 1/(1+4 * a);
    float tmp[SIZE];
    __m256 _top, _bottom, _left, _right, _resA, _resB, _a, _a_div, _center;
    _a = _mm256_set1_ps(a);
    _a_div = _mm256_set1_ps(a_div);
    for(int k = 0; k < 20; k++){
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j+=8){
                _top = _mm256_loadu_ps(&x[IX(i-1, j)]);
                _bottom = _mm256_loadu_ps(&x[IX(i+1, j)]);
                _left = _mm256_loadu_ps(&x[IX(i, j-1)]);
                _right = _mm256_loadu_ps(&x[IX(i, j+1)]);
                _center = _mm256_loadu_ps(&x0[IX(i, j)]);
                _resA = _mm256_add_ps(_top,_bottom);
                _resB = _mm256_add_ps(_left,_right);
                _resA = _mm256_add_ps(_resA,_resB);
                _resA = _mm256_fmadd_ps(_resA, _a , _center);
                //_resA = _mm256_mul_ps(_resA, _a);
                //_resA = _mm256_add_ps(_resA,_center);
                _resA = _mm256_mul_ps(_resA, _a_div);


                _mm256_storeu_ps(&tmp[IX(i,j)], _resA);
            }
        }
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j++){
                x[IX(i,j)] = tmp[IX(i,j)];
            }
        }
        set_bnd(n, b, x);
    }
}
//#pragma optimize("", on)

void advect(int n, int b, float* d, float* d0, float* u, float* v, float dt){
    float dt0 = dt*n;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            float x = i-(dt0 * u[IX(i,j)]);
            float y = j-(dt0 * v[IX(i,j)]);

            x = MAX(x,0.5f);
            x = MIN(x,n+0.5f);
            y = MAX(y,0.5f);
            y = MIN(y,n+0.5f);

            int i0 = (int)x;
            int j0 = (int)y;
            int i1 = i0 + 1;
            int j1 = j0 + 1;

            float s1 = x-i0;
            float t1 = y - j0;
            float s0 = 1 - s1;
            float t0 = 1 - t1;

            float a = t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)];
            float b =  t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)];

            d[IX(i,j)] = s0 * a + s1 * b;
        }
    }
    set_bnd(n, b, d);
}


void project(int n, float* u, float* v, float* p, float* div){
    __m256 _top, _bottom, _left, _right, _resA, _resB, _quarter, _h,_halfn, _center, _center2, _zero;
    
    float h = -0.5f * 1.0/n;
    float halfn = 0.5f * n;
    _h = _mm256_set1_ps(h);
    _zero = _mm256_set1_ps(0.0f);
    _halfn = _mm256_set1_ps(halfn);
    _quarter = _mm256_set1_ps(0.25f);
    for(int i = 1; i <=n; i++){
        for(int j = 1; j <= n; j+=8){
            //div[IX(i,j)] = h * (u[IX(i+1,j)] - u[IX(i-1,j)] + v[IX(i,j+1)] - v[IX(i,j-1)]);
            //p[IX(i,j  )] = 0;
            _top = _mm256_loadu_ps(&u[IX(i-1, j)]);
            _bottom = _mm256_loadu_ps(&u[IX(i+1, j)]);
            _left = _mm256_loadu_ps(&v[IX(i, j-1)]);
            _right = _mm256_loadu_ps(&v[IX(i, j+1)]);
            _resA = _mm256_sub_ps(_bottom,_top);
            _resB = _mm256_sub_ps(_right,_left);
            _resA = _mm256_add_ps(_resA,_resB);
            _resA = _mm256_mul_ps(_resA,_h);
            
            _mm256_storeu_ps(&div[IX(i,j)], _resA);
            _mm256_storeu_ps(&p[IX(i,j)], _zero);
        }
    }
    set_bnd(n, 0, div);
    set_bnd(n, 0, p);
    
    float tmp[SIZE] = {0.0f};
    for(int k = 0; k<20; k++){
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <= n; j+=8){
             //   tmp[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)]) * 0.25f;
                _top = _mm256_loadu_ps(&p[IX(i-1, j)]);
                _bottom = _mm256_loadu_ps(&p[IX(i+1, j)]);
                _left = _mm256_loadu_ps(&p[IX(i, j-1)]);
                _right = _mm256_loadu_ps(&p[IX(i, j+1)]);
                _center = _mm256_loadu_ps(&div[IX(i,j)]);
                _resA = _mm256_add_ps(_bottom,_top);
                _resB = _mm256_add_ps(_right,_left);
                _resA = _mm256_add_ps(_resA,_resB);
                _resA = _mm256_add_ps(_resA,_center);
                _resA = _mm256_mul_ps(_resA,_quarter);
                _mm256_storeu_ps(&tmp[IX(i,j)], _resA);
            }
        }
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j++){
                p[IX(i,j)] = tmp[IX(i,j)];
            }
        }
        set_bnd(n,0,p);
    }

    for(int i = 1; i <=n; i++){
        for(int j = 1; j <= n; j+=8){
            //BUG
            //u[IX(i,j)] -= halfn * (p[IX(i+1, j)] - p[IX(i-1,j)]);
            //v[IX(i,j)] -= halfn * (p[IX(i, j+1)] - p[IX(i,j-1)]);
            //*
            _top = _mm256_loadu_ps(&p[IX(i-1, j)]);
            _bottom = _mm256_loadu_ps(&p[IX(i+1, j)]);
            _center = _mm256_loadu_ps(&u[IX(i,j)]);
            _resA = _mm256_sub_ps(_bottom,_top);
            _resA = _mm256_mul_ps(_resA,_halfn);
            _resA = _mm256_sub_ps(_center,_resA);
            _mm256_storeu_ps(&u[IX(i,j)], _resA);


            _left = _mm256_loadu_ps(&p[IX(i, j-1)]);
            _right = _mm256_loadu_ps(&p[IX(i, j+1)]);
            _center2 = _mm256_loadu_ps(&v[IX(i,j)]);
            _resB = _mm256_sub_ps(_right,_left);
            _resB = _mm256_mul_ps(_resB,_halfn);
            _resB = _mm256_sub_ps(_center2,_resB);
            _mm256_storeu_ps(&v[IX(i,j)], _resB);


        }
    }
    set_bnd(n, 1, u);
    set_bnd(n,2,v);
}
#endif
