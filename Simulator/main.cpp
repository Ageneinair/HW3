//
//  main.cpp
//  Simulator
//
//  Created by CHENZHANG WANG on 11/19/18.
//  Copyright © 2018 CHENZHANG WANG. All rights reserved.
//

#include <iostream>
#include <math.h>

double array_minus(double array1[], double array2[]);
double add(double array1[], double array2[]);
double normalize(double array[]);


class cube
{
public:
    double initial_position[27][3]
    double position[27][3];
    double velocity[27][3];
    double acceration[27][3];
    double force[27][3];
}

int main(int argv, char *argv)
{
    cube  c;
    
    for(int x=0; x<3; x++)
    {
        for(int y=0; y<3; y++)
        {
            for(int z=0; z<3; z++)
            {
                c.position[x+y+z][0] = x;
                c.position[x+y+z][1] = y;
                c.position[x+y+z][2] = z;
            }
        }
    }
    c.initial_position = c.position
}

double array_minus(double array1[], double array2[])
{
    length = sizeof(array);
    double difference[length];
    for(int i = 0; i < length; i++)
    {
        difference[i] = array1[i] - array[i];
    }
    return difference;
}

double add(double array1[], double array2[])
{
    length = sizeof(array);
    double sum[length];
    for(int i = 0; i < length; i++)
    {
        sum[i] = array1[i] + array[i];
    }
    return sum;
}

double normalize(double array[])
{
    length = sizeof(array);
    double new_array[length];
    for(int i = 0; i < length; i++)
    {
        sum += pow(array[i],2);
    }
    sum = sqrt(sum);
    for(int i = 0; i < length; i++)
    {
        new_array[i] = array[i] / sum;
    }
    return new_array;
}
