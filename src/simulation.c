#ifndef SIM_ORIGINAL_C
#define SIM_ORIGINAL_C 

#include "immintrin.h"
//#TODO: how to define this in main
#define N 200
#define SIZE (N + 2) * (N + 2)
#define IX(i,j) ((j) + (N+2)*(i))
#define SWAP(x0,x) {float *tmp = x0; x0=x;x=tmp;}
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

void add_source(int n, float* x, float* s, float dt){
    int size = (n+2) * (n+2);
    for(int i = 0; i < size; i++){
        x[i] += dt * s[i];  
    }
}

void set_bnd(int n, int b, float* x){
    for(int i = 1; i < n; i++){
        x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
        x[IX(n+1,i)] = b==1 ? -x[IX(n,i)] : x[IX(n,i)];
        x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,n+1)] = b==2 ? -x[IX(i,n)] : x[IX(i,n)];
    }
    x[IX(0  ,0  )] = 0.5 * (x[IX(1,0  )] + x[IX(0  ,1)]);
    x[IX(0  ,n+1)] = 0.5 * (x[IX(1,n+1)] + x[IX(0  ,n)]);
    x[IX(n+1,0  )] = 0.5 * (x[IX(n,0  )] + x[IX(n+1,1)]);
    x[IX(n+1,n+1)] = 0.5 * (x[IX(n,n+1)] + x[IX(n+1,n)]);
}

void diffuse ( int n, int b, float * x, float * x0, float diff, float dt )
{
    int i, j, k;
    float a=dt*diff*n*n;
    for ( k=0 ; k<20 ; k++ ) {
        for ( i=1 ; i<=n ; i++ ) {
            for ( j=1 ; j<=n ; j++ ) {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
            }
        }
        set_bnd ( n, b, x );
    }
}

void diffuse_pipelined(int n, int b, float* x, float* x0, float diff, float dt){
    float a = dt * diff * n * n;
    float a_div = 1/(1+4 * a);
    float tmp[SIZE];
    for(int k = 0; k < 20; k++){
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j+=4){
                float sum1= x[IX(i - 1,j  )] + x[IX(i+1,j  )] + x[IX(i,j-1)] + x[IX(i,j+1)];
                float sum2= x[IX(i - 1,j+1)] + x[IX(i+1,j+1)] + x[IX(i,j  )] + x[IX(i,j+2)];
                float sum3= x[IX(i - 1,j+2)] + x[IX(i+1,j+2)] + x[IX(i,j+1)] + x[IX(i,j+3)];
                float sum4= x[IX(i - 1,j+3)] + x[IX(i+1,j+3)] + x[IX(i,j+2)] + x[IX(i,j+4)];

                tmp[IX(i,j  )] = (x0[IX(i,j  )] + (a * sum1)) * a_div;
                tmp[IX(i,j+1)] = (x0[IX(i,j+1)] + (a * sum2)) * a_div;
                tmp[IX(i,j+2)] = (x0[IX(i,j+2)] + (a * sum3)) * a_div;
                tmp[IX(i,j+3)] = (x0[IX(i,j+3)] + (a * sum4)) * a_div;
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

void advect(int n, int b, float* d, float* d0, float* u, float* v, float dt){
    float dt0 = dt*n;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            float x = i-(dt0 * u[IX(i,j)]);
            float y = j-(dt0 * v[IX(i,j)]);

            if(x < 0.5f){ x = 0.5;}
            if(x > n + 0.5f){x = n + 0.5;}
            int i0 = (int)x;
            int i1 = i0 + 1;

            if(y<0.5){y = 0.5;}
            if(y>n+0.5){y = n+0.5;}
            int j0 = (int)y;
            int j1 = j0 + 1;

            float s1 = x-i0;
            float s0 = 1 - s1;

            float t1 = y - j0;
            float t0 = 1 - t1;

            d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) + s1 * (t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]);
            d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) + s1 * (t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]); 
        }
    }
    set_bnd(n, b, d);
}

void project( int n, float * u, float * v, float * p, float * div ){
    int i, j, k;
    float h;
    h = 1.0/n;
    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+
            v[IX(i,j+1)]-v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }
            set_bnd ( n, 0, div ); set_bnd ( n, 0, p );
    for ( k=0 ; k<20 ; k++ ) {
        for ( i=1 ; i<=n ; i++ ) {
            for ( j=1 ; j<=n ; j++ ) {
                p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+p[IX(i,j-1)]+p[IX(i,j+1)])/4;
            }
        }
        set_bnd ( n, 0, p );
    }

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    set_bnd ( n, 1, u ); set_bnd ( n, 2, v );
}

void project_pipelined(int n, float* u, float* v, float* p, float* div){
    float h = -0.5f * 1.0/n;
    for(int i = 1; i <=n; i++){
        for(int j = 1; j <= n; j++){
            div[IX(i,j)] = h * (u[IX(i+1,j)] - u[IX(i-1,j)] + v[IX(i,j+1)] - v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }
    set_bnd(n, 0, div);
    set_bnd(n, 0, p);
    
    float tmp[SIZE] = {0.0f};
    for(int k = 0; k<20; k++){
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <= n; j++){
            //p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)]) * 0.25f;
            tmp[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)]) * 0.25f;
            }
        }
        for(int i = 1; i <=n; i++){
            for(int j = 1; j <=n; j++){
                p[IX(i,j)] = tmp[IX(i,j)];
            }
        }
        set_bnd(n,0,p);
    }

    float halfn = 0.5f * n;
    for(int i = 1; i <=n; i++){
        for(int j = 1; j <= n; j++){
            u[IX(i,j)] -= halfn * (p[IX(i+1, j)] - p[IX(i-1,j)]);
            v[IX(i,j)] -= halfn * (p[IX(i, j+1)] - p[IX(i,j-1)]);
        }
    }
    set_bnd(n, 1, u);
    set_bnd(n,2,v);
}
#endif
