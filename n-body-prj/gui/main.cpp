#include <windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <math.h>

#include "GSimulation.hpp"

HDC hDC;
HGLRC hRC;
HWND hWnd;
HINSTANCE hInstance;

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
HANDLE hWorkerThread = NULL;

struct ThreadData {
    int body_num;
    int time;
    real_type area_lim;
};

void EnableOpenGL(HWND hWnd, HDC* hDC, HGLRC* hRC) {
    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR), 1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
        PFD_TYPE_RGBA, 24, 0,0,0,0,0,0,
        0,0,0,0,0,0,0,
        32, 0, 0,
        PFD_MAIN_PLANE, 0, 0,0,0
    };

    *hDC = GetDC(hWnd);
    int iFormat = ChoosePixelFormat(*hDC, &pfd);
    SetPixelFormat(*hDC, iFormat, &pfd);
    *hRC = wglCreateContext(*hDC);
    wglMakeCurrent(*hDC, *hRC);
}

void DisableOpenGL(HWND hWnd, HDC hDC, HGLRC hRC) {
    wglMakeCurrent(NULL, NULL);
    wglDeleteContext(hRC);
    ReleaseDC(hWnd, hDC);
}

struct Sphere {
    GLfloat ambient[4];
    GLfloat diffuse[4];
    GLfloat specular[4];
    GLfloat shininess;
};

DWORD WINAPI WorkerThread(LPVOID lpParam) {
    ThreadData* data = (ThreadData*)lpParam;
    GSimulation sim;
    sim.set_number_of_particles(data->body_num);
    sim.set_sim_time(data->time);
    sim.set_area_size(data->area_lim);

    EnableOpenGL(hWnd, &hDC, &hRC);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.0, 600.0 / 600.0, 1.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, 600, 600);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    GLfloat light_pos[] = { 0.0f, 5.0f, 5.0f, 1.0f };
    GLfloat light_amb[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat light_dif[] = { 0.9f, 0.9f, 0.9f, 1.0f };
    GLfloat light_spe[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_amb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_dif);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_spe);

    Sphere spheres[3] = {
        { {0.2f, 0.1f, 0.0f, 1.0f}, {1.0f, 0.5f, 0.2f, 1.0f}, {1.0f, 1.0f, 1.0f, 1.0f}, 30.0f },
        { {0.0f, 0.2f, 0.2f, 1.0f}, {0.2f, 1.0f, 1.0f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, 60.0f },
        { {0.1f, 0.0f, 0.2f, 1.0f}, {0.7f, 0.2f, 1.0f, 1.0f}, {1.0f, 0.6f, 1.0f, 1.0f}, 90.0f }
    };

    sim.start([&]{
            size_t parts_num = sim.get_part_num();
            Particle* parts = sim.get_parts();
            real_type area_lim = sim.get_area_lim();

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glLoadIdentity();

            glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

            GLUquadric* quad = gluNewQuadric();

            for(size_t i = 0; i < parts_num; ++i) {
                if (parts[i].pos[0] < area_lim && parts[i].pos[1] < area_lim) {
                        glPushMatrix();
                        glTranslatef(((parts[i].pos[0] / area_lim) - 0.5) * 10, ((parts[i].pos[1] / area_lim) - 0.5) * 10, ((parts[i].pos[2] / area_lim) - 0.5 * 10));
                        glMaterialfv(GL_FRONT, GL_AMBIENT,   spheres[i % 3].ambient);
                        glMaterialfv(GL_FRONT, GL_DIFFUSE,   spheres[i % 3].diffuse);
                        glMaterialfv(GL_FRONT, GL_SPECULAR,  spheres[i % 3].specular);
                        glMaterialf (GL_FRONT, GL_SHININESS, spheres[i % 3].shininess);

                        gluSphere(quad, 0.005 * parts[i].mass, 40, 40);
                        glPopMatrix();

                    }
            }
            gluDeleteQuadric(quad);
            SwapBuffers(hDC);
    });

    DisableOpenGL(hWnd, hDC, hRC);
    return 0;
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int) {
    hInstance = hInst;

    WNDCLASS wc = { 0 };
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = "GLSphereClass";
    RegisterClass(&wc);
    hWnd = CreateWindow("GLSphereClass", "Multiple Unique Spheres",
    WS_OVERLAPPEDWINDOW, 100, 100, 600, 600,
    NULL, NULL, hInstance, NULL);

    ShowWindow(hWnd, SW_SHOW);

    MSG msg;
    BOOL bQuit = FALSE;
    float t = 0.0f;

    int argc;
    wchar_t **argv = CommandLineToArgvW(GetCommandLineW(), &argc);

    if (argv == NULL) {
        std::cerr << "Error parsing command line" << std::endl;
        return 1;
    }

    int pointCount = 1000;
    int time = 10;
    real_type area = 0.1;
    wchar_t* endPtr;
    if (argc != 4) {
        std::cout << "Use: n_body_simulation.exe PARTICLES_NUM SIM_TIME AREA_LIM\n";
        return 0;
    }
    pointCount = std::wcstol(argv[1], &endPtr, 10);
    time = std::wcstol(argv[2], &endPtr, 10);
    area = std::wcstod(argv[3], &endPtr);


    ThreadData data{pointCount, time, area};

    hWorkerThread = CreateThread(NULL, 0, WorkerThread, &data, 0, NULL);
    if (hWorkerThread == NULL) {
        MessageBox(NULL, "Failed to create worker thread", "Error!", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    while(GetMessage(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }


    DestroyWindow(hWnd);
    return msg.wParam;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
switch (message) {
case WM_CLOSE:
    PostQuitMessage(0);
    return 0;
case WM_DESTROY:
    return 0;
default:
    return DefWindowProc(hWnd, message, wParam, lParam);
}
}