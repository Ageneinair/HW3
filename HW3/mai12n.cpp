//
//  main.cpp
//  HW3
//
//  Created by CHENZHANG WANG on 11/19/18.
//  Copyright © 2018 CHENZHANG WANG. All rights reserved.
//


#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;

double norm(double array[]);

class cube
{
public:
    double initial_position[8][3] = { 0 };
    double position[8][3] = { 0 };
    double velocity[8][3] = { 0 };
    double acceration[8][3] = { 0 };
    double force[8][3] = { 0 };
    double link[28][3] = { 0 };
};

double get_fitness()
{
    const double g = 9.8;
    const double m = 0.1;
    const int spring_constant = 10000;
    const double l0 = 0.1;
    const double dt = 0.0001;
    const double damping = 1;
    const double initial_height = 0.5;
    const int T = 100000;
    double fitness = 0.0;
    const int n = 2;
    // Initialize a cube
    cube  c;
    int i = 0;
    for (double x=0; x<l0*n; x+=l0)
    {
        for (double y=0; y<l0*n; y+=l0)
        {
            for (double z=0; z<l0*n; z+=l0)
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
    double l_0 = 0;
    double delta_l = 0;
    double direction[3] = { 0 };
    double fl[3] = { 0 };
    
    // Store links
    double parameters[28][3]; //a b c
    double links[n*n*n][n*n*n] = {0}; //parameters index
    int links_index = 0;
    for(int i=0; i<n*n*n; i++)
    {
        for(int j=i+1; j<n*n*n; j++)
        {
            links[i][j] = links_index;
            links_index += 1;
        }
    }
    
    //double ep_ground[T] = {0};
    //double ep_spring[T] = {0};
    //double ek[T] = {0};
    //double ep_gravity[T] = {0};
    
    while (1)
    {
        cout<<c.position[0][2]<<endl;
        // Set end point
        if (t >= T)
            break;
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
                    if(l_0 < 2*l0)
                    {
                        //cout << delta_initial_p[0] << endl;
                        int link_index = links[i][j];
                        l_0 = (parameters[link_index][0] + parameters[link_index][1]*sin(t+parameters[link_index][2]));
                        delta_l = norm(delta_p) - l_0;
                        for (int k = 0; k < 3; k++)
                        {
                            direction[k] = delta_p[k] / norm(delta_p);
                            fl[k] = -1 * spring_constant * delta_l * direction[k];
                            c.force[i][k] += fl[k];
                        }
                        //ep_spring[T] += pow(delta_l,2) * k / 2;
                    }
                }
            }
            // Ground force
            if (c.position[i][2] < 0)
            {
                // Vertical
                fz = 10000 * pow(c.position[i][2], 2);
                //ep_ground[t] = pow(c.position[i][2],2) * k / 2;
                //Horizontal
                fxy = sqrt(pow(c.force[i][0], 2) + pow(c.force[i][1], 2));
                vxy = sqrt(pow(c.velocity[i][0], 2) + pow(c.velocity[i][1], 2));
                // Static friction
                if (vxy == 0)
                {
                    c.force[i][2] += fz;
                }
                // Dynamic friction
                else
                {
                    fx = -1 * fz * c.velocity[i][0] / vxy;
                    fy = -1 * fz * c.velocity[i][1] / vxy;
                    c.force[i][2] += fx;
                    c.force[i][2] += fy;
                    c.force[i][2] += fz;
                }
            }
            
            //cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
        }
        //cout << c.force[0][2] << endl;
        
        // Update parameters
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                c.acceration[i][j] = c.force[i][j] / m;
                c.velocity[i][j] += c.acceration[i][j] * dt;
                c.velocity[i][j] *= damping;
                c.position[i][j] += c.velocity[i][j] * dt;
                
            }
            
        }
        t += 1;
        //cout << t << endl;
        
    }
    //cout << c.position[0][1] << endl;
    double final_position[2] = { 0 };
    double first_position[2] = { 0 };
    for (int i = 0; i < 8; i++)
    {
        final_position[0] += c.position[i][0]/8;
        final_position[1] += c.position[i][1]/8;
        first_position[0] += c.initial_position[i][0] / 8;
        first_position[1] += c.initial_position[i][1] / 8;
    }
    fitness = sqrt(pow(final_position[0] - first_position[0],2)+ pow(final_position[1] - first_position[1], 2));
    return fitness;
}

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

int main()
{
    get_fitness();
    //cout << get_fitness() << endl;
    cin.get();
    return 0;
}
