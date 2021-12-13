#pragma once
/*
* Low storage Runge-kutta
* M.H. Carpenter and C.A. Kennedy. 
* Fourth-order 2N-storage Runge-Kutta schemes. 
* TechnicalReport NASA TM-109112,
* 4th order 5-4 solution 3
*
* To integrate in time from time `t` to `t + dt`, compute:
*
* for (int k = 0; k < 5; ++k) {
*     // Set rates
*     dy = a[k]*dy + g(y,t+c[k]*dt); 
* 
*     // Update fields
*     y = y + dt*b[k]*dy;
* }
*
%  g: g(y, t) function to integrate in time
%  y: Solution vector at current time t
% dt: Time step
*
*/ 

template <typename rk4Functor, typename callbackFunctor,
          typename Tv=double>
void lsrk4(rk4Functor& rk4, callbackFunctor& callback,
           const Tv T, const Tv dt) {

                Tv rk4_a[5] = {0.0,
                                -567301805773.0/1357537059087,
                                -2404267990393.0/2016746695238,
                                -3550918686646.0/2091501179385,
                                -1275806237668.0/842570457699};
                Tv rk4_b[5] = {1432997174477.0/9575080441755,
                                5161836677717.0/13612068292357,
                                1720146321549.0/2090206949498,
                                3134564353537.0/4481467310338,
                                2277821191437.0/14882151754819};
                Tv rk4_c[5] = {0.0,
                                             1432997174477.0/9575080441755,
                                             2526269341429.0/6820363962896,
                                             2006345519317.0/3224310063776,
                                             2802321613138.0/2924317926251};

        callback(0, 0, dt);
        int num_steps = (int) ((T - 0.5 * dt) / dt);
        for (int step = 0; step < num_steps; ++step) {
                Tv t = step * dt;
                for (int k = 0; k < 5; ++k) {
                        rk4.rates(t + rk4_c[k] * dt, rk4_a[k]);
                        rk4.update(dt * rk4_b[k]);
                }
                callback(step + 1, t + dt, dt);

        }

        // Take final correction step to reach the desired final time
        Tv new_dt = T - (num_steps) * dt;
        Tv t = num_steps * dt;
        if (new_dt > 0) {
                for (int k = 0; k < 5; ++k) {
                                rk4.rates(t + rk4_c[k] * new_dt, rk4_a[k]);
                                rk4.update(new_dt * rk4_b[k]);
                }
        }
        callback(num_steps, T, dt, 1);
}
