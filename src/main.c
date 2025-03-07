#include <bits/time.h>
#include <stdint.h>
#include <unistd.h>
#include "raylib.h"
#include <time.h>
#include "simulation_vector.c"


//TODO:
//Reset Button
//UI to modify
//None Square size
//Faster Render
//Pause/Rewind/Step
//Zoom | camera or some way to only show pieces of it, not the entire thing
//Adjustable sensitivity
//Walls/Objects

#define N 200
#define SIZE (N + 2) * (N + 2)
#define IX(i,j) ((j) + (N+2)*(i))
#define SWAP(x0,x) {float *tmp = x0; x0=x;x=tmp;}
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


//TODO: memory arena maybe? just to say I've done it

//TODO: otherwise preset size in struct so it's one block of memory;
typedef struct perfLog{
     int *times;
     int *ops;
     //idx and clocks dont really belong here, but whatever for now
     struct timespec startClock;
     struct timespec endClock;
     int idx;

} perfLog;

typedef struct fluidField{
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
} fluidField;


typedef enum{
    DYE,
    FLOW,
}renderMode;

int screenWidth = 1000;
int screenHeight = 1000;
RenderTexture2D texture;
perfLog perf;
renderMode rendering = DYE;

void setup(){
    InitWindow(screenWidth,screenHeight,"Jam Game");
    texture = LoadRenderTexture(screenWidth, screenHeight);
    //TODO: this is just kinda temporary and jank
    perf.idx = 0;
    perf.times = malloc(10 * sizeof(int));
    perf.ops = malloc(10 * sizeof(long));
    perf.ops[0] = 20 * N * N * 5;
    perf.ops[1] = 20 * N * N;
    perf.ops[2] = (4 * N * N) + (20 * N * N * 5)  + (N * N * 6);
    perf.ops[3] = (5 * perf.ops[0]) + (5 * perf.ops[1]) + perf.ops[2];
    perf.ops[4] = N * N;
    //SetTargetFPS(60);
}


void reset(fluidField field){
    for(int i = 0; i < SIZE; i++){
        field.u[i] = 0.0f;
        field.v[i] = 0.0f;
        field.dens_red[i] = 0.0f;
        field.dens_green[i] = 0.0f;
        field.dens_blue[i] = 0.0f;
        field.u_prev[i] = 0.0f;
        field.v_prev[i] = 0.0f;
        field.dens_red_prev[i] = 0.0f;
        field.dens_green_prev[i] = 0.0f;
        field.dens_blue_prev[i] = 0.0f;
    }
}


void addDensity(int i, int j, float*x,  float amount){
    x[IX(i,j)] += amount;
}
void addVelocity(int i, int j, float*x,float*y,  float amountX, float amountY){
    x[IX(i,j)] += amountX;
    y[IX(i,j)] += amountY;
}


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
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.startClock);
    diffuse(n, 1, u, u0, visc, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.endClock);
    perf.times[0] = perf.endClock.tv_nsec - perf.startClock.tv_nsec;
    SWAP(v0,v);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.startClock);
    diffuse(n, 2 ,v, v0, visc, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.endClock);
    perf.times[2] = perf.endClock.tv_nsec - perf.startClock.tv_nsec;
    project(n, u, v, u0, v0);
    SWAP(u0,u);
    SWAP(v0,v);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.startClock);
    advect(n, 1, u, u0, u0, v0, dt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.endClock);
    perf.times[1] = perf.endClock.tv_nsec - perf.startClock.tv_nsec;
    advect(n, 2, v, v0, u0, v0, dt);
    project(n, u, v, u0, v0);
}




void DrawScene(float* red, float* green, float* blue){
    BeginDrawing();
    {
        float sensitivityMultiplier = 1000.0f;
        int max = 255;
        int min = 0;
        int offset = 0;
        if(rendering == FLOW){
            sensitivityMultiplier = 1000.0f;
            max = 127;
            min = -128;
            offset = 128;
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.startClock);
        int scale = screenHeight / N;
        for(int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
                char redVal = MAX(min,MIN(max, (int)(red[IX(i,j)] * sensitivityMultiplier))) + offset;
                char greenVal = MAX(min,MIN(max, (int)(green[IX(i,j)] *sensitivityMultiplier)));// + offset;
                char blueVal = MAX(min,MIN(max, (int)(blue[IX(i,j)] * sensitivityMultiplier))) + offset;
                Color denityColor = {.r= redVal, .g = greenVal, .b = blueVal, .a = 255};
                DrawRectangle((i - 1) * scale, (j-1) * scale, scale, scale, denityColor);
            }
        }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.endClock);
        perf.times[4] = perf.endClock.tv_nsec - perf.startClock.tv_nsec;

        DrawFPS(10,10);
        DrawText(TextFormat("time: %2.06f", perf.times[perf.idx] / (float)1000000000), 10, 30, 20, GREEN);
        long clockSpeed = 3800000000;
        float clockCycles = (perf.times[perf.idx] / (float)1000000000) * clockSpeed;
        DrawText(TextFormat("MCycles: %6.03f", clockCycles / 1000000), 10, 50, 20, GREEN);
        DrawText(TextFormat("Operations/Cycle: %02.01f", perf.ops[perf.idx] / clockCycles), 10, 70, 20, GREEN);
        EndTextureMode();


    }
    EndDrawing();
}



//TODO: draw scene only on a delta time?
void GameLoop(){
 
    fluidField field;
    field.u = malloc(SIZE *sizeof(float));
    field.v = malloc(SIZE *sizeof(float));
    field.dens_red = malloc(SIZE *sizeof(float));
    field.dens_green = malloc(SIZE *sizeof(float));
    field.dens_blue = malloc(SIZE *sizeof(float));
    field.u_prev = malloc(SIZE *sizeof(float));
    field.v_prev = malloc(SIZE *sizeof(float));
    field.dens_red_prev = malloc(SIZE *sizeof(float));
    field.dens_green_prev = malloc(SIZE *sizeof(float));
    field.dens_blue_prev = malloc(SIZE *sizeof(float));

    float dt = 1/60.0f;//0.1f;
    float visc = 0.0000148f;
    float diff = 0.0000001f;


    while(!WindowShouldClose()){
        if(IsKeyDown(KEY_ONE)){rendering = DYE;}
        if(IsKeyDown(KEY_TWO)){rendering = FLOW;}
        if(IsKeyDown(KEY_KP_0)){perf.idx = 0;}
        if(IsKeyDown(KEY_KP_1)){perf.idx = 1;}
        if(IsKeyDown(KEY_KP_2)){perf.idx = 2;}
        if(IsKeyDown(KEY_KP_3)){perf.idx = 3;}
        if(IsKeyDown(KEY_KP_4)){perf.idx = 4;}

        for(int i = 0; i < 5;i++){
            addVelocity(N/2 - 60, N/2  + i, field.u_prev, field.v_prev, 30.f,0.0f);
            addVelocity(N/2 + 60, N/2  + i, field.u_prev, field.v_prev, -30.0f,0.0f);
            addDensity(N/2 - 50,  N/2 + i, field.dens_red_prev, 1.0f * N * N / 20 / 800);
            addDensity(N/2 + 50, N/2 + i, field.dens_blue_prev, 1.0f * N * N  / 20 / 800);
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.startClock);
        int startTime = perf.startClock.tv_nsec;
        vel_step(N, field.u, field.v, field.u_prev, field.v_prev, visc, dt);
        dens_step(N, field.dens_red, field.dens_red_prev, field.u, field.v, diff, dt);
        dens_step(N, field.dens_green, field.dens_green_prev, field.u, field.v, diff, dt);
        dens_step(N, field.dens_blue, field.dens_blue_prev, field.u, field.v, diff, dt);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &perf.endClock);
        int endTime = perf.endClock.tv_nsec;
        perf.times[3] = endTime - startTime;
        for(int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
                field.u_prev[IX(i,j)] = 0;
                field.v_prev[IX(i,j)] = 0;
                field.dens_red_prev[IX(i,j)] = 0;
                field.dens_green_prev[IX(i,j)] = 0;
                field.dens_blue_prev[IX(i,j)] = 0;
            }
        }
        if(rendering == DYE){
            DrawScene(field.dens_red, field.dens_green, field.dens_blue);
        }else{
            float empty[SIZE] = {-128.0f};
            DrawScene(field.u, empty, field.v);
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
