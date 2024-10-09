#include "udf.h"

// Define a generalized velocity profile with adjustable parameters
DEFINE_PROFILE(velocity_profile, thread, index)
{
    float x[ND_ND];      // Coordinates array to hold spatial positions
    face_t f;            // Face index for looping through each face in the thread
    real t = CURRENT_TIME;  // Current time in the simulation

    // Begin loop over each face in the specified thread
    begin_f_loop(f, thread)
    {
        F_CENTROID(x, f, thread);  // Get the centroid coordinates of the face
        
        // Define time intervals and polynomial expressions for velocity profile
        if (t > 0 && t <= 1)
        {
            F_PROFILE(f, thread, index) = (2 - 766805 * (x[0] * x[0]) - 766805 * (x[2] * x[2])) * (
                3680.169708 * pow(t, 13) + 63303.80158 * pow(t, 12) - 488379.7891 * pow(t, 11) +
                1452464.048 * pow(t, 10) - 2430160.03 * pow(t, 9) + 2575471.32 * pow(t, 8) -
                1814743.614 * pow(t, 7) + 862641.3857 * pow(t, 6) - 273760.3681 * pow(t, 5) +
                55808.39602 * pow(t, 4) - 6708.998621 * t * t * t + 388.7979179 * t * t -
                5.131067551 * t + 0.3892026544);
        }
        else if (t > 1 && t <= 2)
        {
            real time_shift = t - 1;
            F_PROFILE(f, thread, index) = (2 - 766805 * (x[0] * x[0]) - 766805 * (x[2] * x[2])) * (
                3680.169708 * pow(time_shift, 13) + 63303.80158 * pow(time_shift, 12) - 
                488379.7891 * pow(time_shift, 11) + 1452464.048 * pow(time_shift, 10) -
                2430160.03 * pow(time_shift, 9) + 2575471.32 * pow(time_shift, 8) -
                1814743.614 * pow(time_shift, 7) + 862641.3857 * pow(time_shift, 6) -
                273760.3681 * pow(time_shift, 5) + 55808.39602 * pow(time_shift, 4) -
                6708.998621 * time_shift * time_shift * time_shift +
                388.7979179 * time_shift * time_shift - 5.131067551 * time_shift + 
                0.3892026544);
        }
        else if (t > 2 && t <= 3)
        {
            real time_shift = t - 2;
            F_PROFILE(f, thread, index) = (2 - 766805 * (x[0] * x[0]) - 766805 * (x[2] * x[2])) * (
                3680.169708 * pow(time_shift, 13) + 63303.80158 * pow(time_shift, 12) - 
                488379.7891 * pow(time_shift, 11) + 1452464.048 * pow(time_shift, 10) -
                2430160.03 * pow(time_shift, 9) + 2575471.32 * pow(time_shift, 8) -
                1814743.614 * pow(time_shift, 7) + 862641.3857 * pow(time_shift, 6) -
                273760.3681 * pow(time_shift, 5) + 55808.39602 * pow(time_shift, 4) -
                6708.998621 * time_shift * time_shift * time_shift +
                388.7979179 * time_shift * time_shift - 5.131067551 * time_shift + 
                0.3892026544);
        }
    }
    end_f_loop(f, thread);  // End loop over all faces
}
