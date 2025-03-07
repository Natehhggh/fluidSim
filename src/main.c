#include <bits/time.h>
#include <stdint.h>
#include <unistd.h>
#include "raylib.h"
#include <time.h>
#include "./simulation_vector.c"


//TODO:
//Reset Button
//UI to modify
//None Square size
//Faster Render
//Pause/Rewind/Step
//Zoom
//Walls/Objects

int screenWidth = 1000;
int screenHeight = 1000;
#define N 200
#define SIZE (N + 2) * (N + 2)
RenderTexture2D texture;

struct timespec startClock, endClock;
int perfTimes[10] = {0};
long perfOps[10] = {0};
int perfId = 0;

/*typedef struct fluidField{
    float *u;
    float *v;
    float *dens_red;
    float *dens_green;
    float *dens_blue;
    float *u_prev;
    float *v_prev;
    float *dens_red_prev;
    float *dens_green_prev;
    float *dens_blue_prev;
} fluidField;*/


#define IX(i,j) ((j) + (N+2)*(i))
#define SWAP(x0,x) {float *tmp = x0; x0=x;x=tmp;}
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
typedef enum{
    DYE,
    FLOW,
}renderMode;

void addDensity(int i, int j, float*x,  float amount){
    x[IX(i,j)] += amount;
}
void addVelocity(int i, int j, float*x,float*y,  float amountX, float amountY){
    x[IX(i,j)] += amountX;
    y[IX(i,j)] += amountY;
}

renderMode rendering = DYE;




void dens_step(int n, float* x, float* x0, float* u, float* v, float diff, float dt){
    add_source(n, x, x0, dt);
    SWAP(x0, x);
    diffuse(n, 0, x, x0, diff, dt);
    SWAP(x0, x);
    advect(n, 0, x, x0, u, v, dt);
}

void vel_step(int n, float* u, float* v, float* u0, float* v0, float visc, float dt){
    add_source(n, u, u0, dt);
    add_source(n,v,v0,dt);
    SWAP(u0,u);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startClock);
    diffuse(n, 1, u, u0, visc, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endClock);
    perfTimes[0] = endClock.tv_nsec - startClock.tv_nsec;
    SWAP(v0,v);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startClock);
    diffuse(n, 2 ,v, v0, visc, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endClock);
    perfTimes[2] = endClock.tv_nsec - startClock.tv_nsec;
    project(n, u, v, u0, v0);
    SWAP(u0,u);
    SWAP(v0,v);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startClock);
    advect(n, 1, u, u0, u0, v0, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endClock);
    perfTimes[1] = endClock.tv_nsec - startClock.tv_nsec;
    advect(n, 2, v, v0, u0, v0, dt);
    project(n, u, v, u0, v0);
}



void setup(){
    InitWindow(screenWidth,screenHeight,"Jam Game");
    texture = LoadRenderTexture(screenWidth, screenHeight);
    //TODO: this is just kinda temporary and jank
    perfOps[0] = 20 * N * N * 5;
    perfOps[1] = 20 * N * N;
    perfOps[2] = (4 * N * N) + (20 * N * N * 5)  + (N * N * 6);
    perfOps[3] = (5 * perfOps[0]) + (5 * perfOps[1]) + perfOps[2];
    //SetTargetFPS(60);
}

void DrawScene(float* red, float* green, float* blue){
    BeginDrawing();
    {
        float sensitivityMultiplier = 1000.0f;
        int max = 255;
        int min = 0;
        int offset = 0;
        int offset_blue = 0;
        if(rendering == FLOW){
            sensitivityMultiplier = 1000.0f;
            max = 127;
            min = -128;
            offset = 128;
        }
        int scale = screenHeight / N;
        for(int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
                char redVal = MAX(min,MIN(max, (int)(red[IX(i,j)] * sensitivityMultiplier))) + offset;
                char greenVal = MAX(min,MIN(max, (int)(green[IX(i,j)] *sensitivityMultiplier))) + offset;
                char blueVal = MAX(min,MIN(max, (int)(blue[IX(i,j)] * sensitivityMultiplier)));// + offset;
                Color denityColor = {.r= redVal, .g = greenVal, .b = blueVal, .a = 255};
                DrawRectangle((i - 1) * scale, (j-1) * scale, scale, scale, denityColor);
            }
        }
        DrawFPS(10,10);
        DrawText(TextFormat("time: %2.06f", perfTimes[perfId] / (float)1000000000), 10, 30, 20, GREEN);
        long clockSpeed = 3800000000;
        float clockCycles = (perfTimes[perfId] / (float)1000000000) * clockSpeed;
        DrawText(TextFormat("MCycles: %6.03f", clockCycles / 1000000), 10, 50, 20, GREEN);
        DrawText(TextFormat("Operations/Cycle: %02.01f", perfOps[perfId] / clockCycles), 10, 70, 20, GREEN);
        EndTextureMode();


    }
    EndDrawing();
}



//TODO: draw scene only on a delta time?
void GameLoop(){
    float *u __attribute__((aligned (64))) = malloc(SIZE *sizeof(float));
    float *v __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_red __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_green __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_blue __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *u_prev __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *v_prev __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_red_prev __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_green_prev __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));
    float *dens_blue_prev __attribute__((aligned (64)))= malloc(SIZE *sizeof(float));

    float dt = 1/60.0f;//0.1f;
    float visc = 0.0000148f;
    float diff = 0.0000000f;
    bool dye = false;


    while(!WindowShouldClose()){
        if(IsKeyDown(KEY_ONE)){rendering = DYE;}
        if(IsKeyDown(KEY_TWO)){rendering = FLOW;}
        if(IsKeyDown(KEY_KP_0)){perfId = 0;}
        if(IsKeyDown(KEY_KP_1)){perfId = 1;}
        if(IsKeyDown(KEY_KP_2)){perfId = 2;}
        if(IsKeyDown(KEY_KP_3)){perfId = 3;}
        if(IsKeyPressed(KEY_KP_9)){dye = !dye;}

        for(int i = 0; i < 5;i++){
            addVelocity(N/2 - 60, N/2  + i, u_prev, v_prev, 50.f,0.0f);
            addVelocity(N/2 + 50, N/2  + i, u_prev, v_prev, -30.0f,0.0f);
            addDensity(N/2 - 50,  N/2 + i, dens_red_prev, 1.0f * N * N / 20 / 500);
            addDensity(N/2 + 50, N/2 + i, dens_blue_prev, 1.0f * N * N  / 20 / 500);
        }

        /*for(int i = 1; i <= N /2;i++){
            addVelocity(5,  i, u_prev, v_prev, 0.1f,0.0f);
            if(i%5 == 0 && dye){
            addDensity(N/2 - 50,  i, dens_red_prev, 1.0f * N * N / 20 / 30000);
            }
        }*/
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startClock);
        int startTime = startClock.tv_nsec;
        vel_step(N, u, v, u_prev, v_prev, visc, dt);
        dens_step(N, dens_red, dens_red_prev, u, v, diff, dt);
        dens_step(N, dens_green, dens_green_prev, u, v, diff, dt);
        dens_step(N, dens_blue, dens_blue_prev, u, v, diff, dt);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endClock);
        int endTime = endClock.tv_nsec;
        perfTimes[3] = endTime - startTime;
        for(int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
                u_prev[IX(i,j)] = 0;
                v_prev[IX(i,j)] = 0;
                dens_red_prev[IX(i,j)] = 0;
                dens_green_prev[IX(i,j)] = 0;
                dens_blue_prev[IX(i,j)] = 0;
            }
        }
        if(rendering == DYE){
            DrawScene(dens_red, dens_green, dens_blue);
        }else{
            float empty[SIZE] = {-128.0f};
            DrawScene(u, v, empty);
        }

    }
}

int main()
{
    setup();
    GameLoop();
    return 0;
}


//#pragma optimize("", on)
