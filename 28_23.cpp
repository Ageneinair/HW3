#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GLTools.h"
#include <GLUT/GLUT.h>
#include <time.h>
#include <iostream>

using namespace std;

class cube
{
public:
    double initial_position[8][3] = { 0 };
    double position[8][3] = { 0 };
    double velocity[8][3] = { 0 };
    double acceration[8][3] = { 0 };
    double force[8][3] = { 0 };
};

double final_result[1000000][8][3] = {0};
static int day = 0;


double norm(double array[])
{
    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        sum += pow(array[i], 2);
    }
    sum = sqrt(sum);
    
    return sum;
}


double CalFrequency() {
    static int count;
    static double save;
    static clock_t last, current;
    double timegap;
    ++count;
    if( count <= 50 )
        return save;
    count = 0;
    last = current;
    current = clock();
    timegap = (current-last)/(double)CLOCKS_PER_SEC;
    save = 50.0/timegap;
    return save;
    
}

void myDisplay(void)
{
    
    double FPS = CalFrequency();
    printf("FPS = %f\n", FPS);
    
    glEnable(GL_DEPTH_TEST);
    // 清除屏幕及深度缓存
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // 深度测试开启，实现遮挡关系
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
    
    //Background color
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    
    //Anti-aliasing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    //    glEnable(GL_MULTISAMPLE);

    //Sight Seting
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, 1, 1, 40000000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(1, 2, 1, 0, 0, 0, 0, 0, 1);
    
    //Groud
    glBegin(GL_QUADS);
    for( int i = -5; i < 5; i++ ) {
        for( int j = -5; j < 5; j++ ) {
            if ((i+j)%2) {
                glColor3f(0.5f, 0.5f, 0.5f);
                glVertex3i(i, j, 0);
                glVertex3i(i+1, j, 0);
                glVertex3i(i+1, j+1, 0);
                glVertex3i(i, j+1, 0);
            }
            else {
                glColor3f(0.1f, 0.1f, 0.1f);
                glVertex3i(i, j, 0);
                glVertex3i(i+1, j, 0);
                glVertex3i(i+1, j+1, 0);
                glVertex3i(i, j+1, 0);

            }
        }
    }
    glEnd();

    //遮挡测试
    glBegin(GL_QUADS);
    glColor3f(1.0f, 0.5f, 0.5f);
    glVertex3i(1, 1,-2);
    glVertex3i(1+1, 1, -2);
    glVertex3i(1+1, 1+1, -2);
    glVertex3i(1, 1+1, -2);
    glEnd();

    //XYZ axis
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f,2.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 2.0f,0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(2.0f, 0.0f, 0.0f);
    glEnd();
    
    //    cout << day << endl;
    //
//        cout << final_result[day][0][0] << ' ' << final_result[day][0][1] << ' '<< final_result[day][0][2] << ' ' << endl;
    


    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 0.0f);
    for ( int i = 0; i < 7; i++){
        for (int j = i+1; j < 8; j ++){
            glVertex3d(final_result[day][i][0],final_result[day][i][1],final_result[day][i][2]);
            glVertex3d(final_result[day][j][0],final_result[day][j][1],final_result[day][j][2]);
        }
    }
    glEnd();
    
    //    //Point Size
    //    glPointSize(5.0f);
    //    glBegin(GL_POINTS);
    //    glVertex3f(0.0f, 0.0f, 0.0f);
    //    glVertex3f(0.0f, 1.0f,0.0f);
    //    glEnd();
    
    glFlush();
    glutSwapBuffers();
}


void myIdle(void) {
    day++;
    while (day % 2 > 0)
        day++;
    if (day>10000)
        day = 0;
    for (int i = 0; i< 100000; i++);
    
    myDisplay();
    
}


int main(int argc, char *argv[])
{
    const double g = 9.8;
    const double m = 0.1;
    const int spring_constant = 10000;
    const int ground_constant =800;
    const double l0 = 0.1;
    const double dt = 0.001;
    const double damping = 1;
    const double initial_height = 0.5;
    const int T = 10000;
    double fitness = 0.0;
    // Initialize a cube
    cube  c;
    int i = 0;
    for (double x = 0; x <= l0; x+=0.1)
    {
        for (double y = 0; y <= l0; y+=0.1)
        {
            for (double z = 0; z <= l0; z+=0.1)
            {
                c.position[i][0] = x;
                c.position[i][1] = y;
                c.position[i][2] = z + initial_height;
                c.initial_position[i][0] = x;
                c.initial_position[i][1] = y;
                c.initial_position[i][2] = z;
                i++;
            }
        }
    }
    
    // Simulator
    int t = 0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double fxy = 0.0;
    double vxy = 0.0;
    double delta_p[3] = { 0 };
    double delta_initial_p[3] = { 0 };
    double d_p[3] = { 0 };
    double l_0 = 0;
    double delta_l = 0;
    double direction[3] = { 0 };
    double fl[3] = { 0 };
    
    //double ep_ground[T] = {0};
    //double ep_spring[T] = {0};
    //double ek[T] = {0};
    //double ep_gravity[T] = {0};
    
    while (t <= T )
    {
//        cout<< c.position[0][2] <<endl;
        // Set end point

        for (int i = 0; i < 8; i++)
        {
            
            // Add gravity
            for (int j = 0; j < 3; j++)
            {
                //cout << i << ' ' << j << endl;
                //cout << m << ' ' << g << endl;
                if (j != 2)
                    c.force[i][j] = 0.0;
                else
                    c.force[i][j] = -1 * m * g;
            }
            
            
            //cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
            
            
            //Sprint force
            for (int j = 0; j < 8; j++)
            {
                if (i != j)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        delta_p[k] = c.position[i][k] - c.position[j][k];
                        delta_initial_p[k] = c.initial_position[i][k] - c.initial_position[j][k];
                    }
                    
                    l_0 = norm(delta_initial_p);
                    //cout << delta_initial_p[0] << endl;
                    delta_l = norm(delta_p) - l_0;
                    
                   // cout << t << ' ' << delta_l << endl;
                    
                    for (int k = 0; k < 3; k++)
                    {
                        direction[k] = delta_p[k] / norm(delta_p);
                        fl[k] = -1 * spring_constant * delta_l * direction[k];
                        c.force[i][k] += fl[k];
                    }
                    //ep_spring[T] += pow(delta_l,2) * k / 2;
                }
            }
            
            // Ground force
            if (c.position[i][2] < 0)
            {
                // Vertical
                fz = -1* ground_constant * c.position[i][2];
                //ep_ground[t] = pow(c.position[i][2],2) * k / 2;
                
                //Horizontal
                fxy = sqrt(pow(c.force[i][0], 2) + pow(c.force[i][1], 2));
                vxy = sqrt(pow(c.velocity[i][0], 2) + pow(c.velocity[i][1], 2));
                //cout << vxy << endl;
                // Static friction
                if (vxy < 0.01 && fxy <= fz * 1)
                {
                    c.force[i][0] = 0;
                    c.force[i][1] = 0;
                    c.force[i][2] += fz;
                }
                // Dynamic friction
                else
                {
                    fx = -0.8 * fz * c.velocity[i][0] / vxy;
                    fy = -0.8 * fz * c.velocity[i][1] / vxy;
                    c.force[i][0] += fx;
                    c.force[i][1] += fy;
                    c.force[i][2] += fz;
                }
            }
            
            //cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
        }
        
        // Update parameters
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                c.acceration[i][j] = c.force[i][j] / m;
                c.velocity[i][j] += c.acceration[i][j] * dt;
                c.velocity[i][j] *= damping;
                c.position[i][j] += c.velocity[i][j] * dt;
                final_result[t][i][j] = c.position[i][j];
            }
        }
        
        t += 1;
        //cout << t << endl;
    }

    
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(400, 400);
    glutCreateWindow("第一个 OpenGL 程序");
    glutDisplayFunc(&myDisplay);
    glutIdleFunc(&myIdle);
    glutMainLoop();
    return 0;
}





