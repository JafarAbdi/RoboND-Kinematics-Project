## Project: Kinematics Pick & Place

[//]: # (Image References)

[DH]: ./misc_images/DH.jpg
[FK]: ./misc_images/FK.png
[P25]: ./misc_images/P25.png
[theta1]: ./misc_images/theta1.png
[theta2]: ./misc_images/theta2.png
[theta3]: ./misc_images/theta3.png

---

### Kinematic Analysis

## 1.DH parameters table

![DH]

Joint i | a(i-1) | alpha(i-1) | d(i) | theta(i)
--- | --- | --- | --- | ---
1 | 0 | 0 | 0.75 | theta1
2 | 0.35 | -90 | 0 | theta2 - 90
3 | 1.25 | 0 | 0 | theta3
4 | -0.054 | -90 | 1.5 | theta4
5 | 0 | 90 | 0 | theta5
6 | 0 | -90 | 0 | theta6
7 | 0 | 0 | 0.453 | 0

* the joint 7 isn't a real joint it's just for transformation from frame 6 to end-effector

I used 'rosrun tf tf_echo' node to calcuate the parameters, for d7 I used the transform between link_5 and left_gripper_finger_link, but without offset in y direction .

to find the homogeneous transformation between two joint I multiplied 4 matrices T(x,a)*R(x,alpha)*T(z,d)*R(z,theta), so to find Ti-1,i; i=1,.....,7 just plug the parameter into this function

```
[                cos(theta(i)),                -sin(theta(i)),               0,               a(i-1)]
[sin(theta(i))*cos(alpha(i-1)), cos(theta(i))*cos(alpha(i-1)),-sin(alpha(i-1)),-sin(alpha(i-1))*d(i)]
[sin(theta(i))*sin(alpha(i-1)), cos(theta(i))*sin(alpha(i-1)), cos(alpha(i-1)), cos(alpha(i-1))*d(i)]
[                            0,                             0,               0,                    1]
```
first I define a dictionary named s and function homo_transformation that take the parameters as input and return the homogeneous transformation and using sympy library

```python
s = {alpha0:     0, a0:     0 , d1:  0.75,
     alpha1: -pi/2, a1:   0.35, d2:     0, q2:q2 - pi/2,
     alpha2:     0, a2:   1.25, d3:     0,
     alpha3: -pi/2, a3: -0.054, d4:   1.5,
     alpha4:  pi/2, a4:      0, d5:     0,
     alpha5: -pi/2, a5:      0, d6:     0,
     alpha6:     0, a6:      0, d7: 0.453, q7: 0}
     
def homo_transformation(alpha, a, d, theta):
    
    return Matrix([[           cos(theta),           -sin(theta),          0,            a],
                   [sin(theta)*cos(alpha), cos(theta)*cos(alpha),-sin(alpha),-sin(alpha)*d],
                   [sin(theta)*sin(alpha), cos(theta)*sin(alpha), cos(alpha), cos(alpha)*d],
                   [                    0,                     0,          0,            1]])
                   
T01 = homo_transformation(alpha0, a0, d1, q1)
T01 = T01.subs(s)
T12 = homo_transformation(alpha1, a1, d2, q2)
T12 = T12.subs(s)
T23 = homo_transformation(alpha2, a2, d3, q3)
T23 = T23.subs(s)
T34 = homo_transformation(alpha3, a3, d4, q4)
T34 = T34.subs(s)
T45 = homo_transformation(alpha4, a4, d5, q5)
T45 = T45.subs(s)
T56 = homo_transformation(alpha5, a5, d6, q6)
T56 = T56.subs(s)
T67 = homo_transformation(alpha6, a6, d7, q7)
T67 = T67.subs(s)
```
after the calculations

T01
```
[cos(q1), -sin(q1), 0,    0] 
[sin(q1),  cos(q1), 0,    0] 
[      0,        0, 1, 0.75] 
[      0,        0, 0,    1]
```
T12
```
[sin(q2),  cos(q2), 0, 0.35]
[      0,        0, 1,    0] 
[cos(q2), -sin(q2), 0,    0] 
[      0,        0, 0,    1]
```
T23
```
[cos(q3), -sin(q3), 0, 1.25] 
[sin(q3),  cos(q3), 0,    0] 
[      0,        0, 1,    0] 
[      0,        0, 0,    1] 
```
T34
```
[ cos(q4), -sin(q4), 0, -0.054] 
[       0,        0, 1,    1.5] 
[-sin(q4), -cos(q4), 0,      0] 
[       0,        0, 0,      1] 
```
T45
```
[cos(q5), -sin(q5), 0, 0] 
[      0,        0,-1, 0] 
[sin(q5),  cos(q5), 0, 0]
[      0,        0, 0, 1] 
```
T56
```
[ cos(q6), -sin(q6), 0, 0]
[       0,        0, 1, 0]
[-sin(q6), -cos(q6), 0, 0]
[       0,        0, 0, 1] 
```
T67
```
[1, 0, 0,     0]
[0, 1, 0,     0]
[0, 0, 1, 0.453]
[0, 0, 0,     1]
```
the homogeneous transformation between base and end-effector can be calculated by multiplying T01*.....T67

T0G = T01 * T12 * T23 * T34 * T45 * T56 * T67; G and 7 are the same and Ti-1,i can be calculated as previously shown and if I want to represent it in urdf represent I should multiply it by correction matrices first rotation about z by pi after that  about y by -pi/2

```python
RC1 = Matrix([[cos(pi),-sin(pi), 0, 0],
              [sin(pi), cos(pi), 0, 0],
              [      0,       0, 1, 0],
              [      0,       0, 0, 1]])
              
RC2 = Matrix([[ cos(-pi/2), 0, sin(-pi/2), 0],
              [          0, 1,          0, 0],
              [-sin(-pi/2), 0, cos(-pi/2), 0],
              [          0, 0,          0, 1]])
RC = simplify(RC1 * RC2)

T07 = simplify(T07 * RC)
```

calculating the forward kinematics after finding the transformation matrices from DH parameters table
![FK]

#### 2. Decouple Inverse Kinematics problem into Inverse Position Kinematics and inverse Orientation Kinematics; doing so derive the equations to calculate all individual joint angles.

Inverse kinematics is the problem of finding the angles that satisfy the given position and orientation of end-effector

## 1. theta 1

I calculated theta1 by first finding the vector P05 which is the vecotr that point from base frame to joint 5

**P05 = P0G - z4*d7** where z4 is the unit vector z of frame 4 and d7 from DH table, P0G is the position which we want our end-effector to in.

after that as shown in the following figure, we calculate theta1 by arctan2(P05_y,P05_x) 

![theta1]

## 2. theta 2

first to calculate theta 2, I found P25 the vector the point from origin of frame 2 to the origin to frame 5 as shown from frame 2 when theta2 is zero 

![P25]

**P25 = R20*(P05 - P02)**; R20 rotation map vector from frame 0 to frame 2 when theta 2 is zero

the long red and green lines are x,y axis respectively of frame 2 when theta2 is zero
![theta2]

theta2 = pi/2 -(beta1 + beta2)

beta1 = arctan2(P25_x,P25_y)

using cosine law

beta2 = arccos((a3^2 + norm(P25)^2 - d^2)/(2*norm(P25)*a2)), where d is calculated using trigonometry d = sqrt(0.054^2 + 1.5^2) and it shown in the previous image, a2 from DH table 

theta2 = pi/2 - (beta1 + beta2)

## 3. theta 3

![theta3]

theta3 = pi/2 - (phi + alpha)

using cosine law 

phi = arccos((a2^2 + d^2 - norm(P25))/(2*a2*d)), d and a2 are as theta2

Links | alpha(i-1) | a(i-1) | d(i-1) | theta(i)
--- | --- | --- | --- | ---
0->1 | 0 | 0 | L1 | qi
1->2 | - pi/2 | L2 | 0 | -pi/2 + q2
2->3 | 0 | 0 | 0 | 0
3->4 |  0 | 0 | 0 | 0
4->5 | 0 | 0 | 0 | 0
5->6 | 0 | 0 | 0 | 0
6->EE | 0 | 0 | 0 | 0


alpha = arctan2(0.054, 0.96)

## 4. theta 4,5,6

first I calculated the general rotation matrix R3G 

-s4*s6 + c4*c5*c6 | -s4*c6 - s6*c4*c5 | -s5*c4 |
-----------------:|------------------:|-------:|
**s5*c6**         |         **-s5*s6**|  **c5**|
**-s4*c5*c6 -s6*c4**|**s4*s6*c5 - c4*c6**|**s4*s5**|

and because I caluclated theta1,2,3 I can find R03, so by comaring I can find the last three angles

theta4 = arctan2(R3G[3,3],-R3G[0,3])

theta5 = arctan2(sqrt(R3G[2,1]^2 + R3G[2,2]^2),R3G[2,3])

theta6 = arctan2(-R3G[2,2],R3G[2,1])

and if sin(theta5) < 0

theta4 = arctan2(-R3G[3,3],R3G[0,3])

theta6 = arctan2(R3G[2,2],-R3G[2,1])

and when theta5 = 0 we have singularity and we can't calculate theta4,6 independenlty we have theta4 + theta6 = const, so one way to solve it(suggested by @gwwang) by letting 

theta4 = 0

theta6 = arctan2(-R36[1,2],-R36[3,2])

### Project Implementation

[Final results](https://www.youtube.com/watch?v=PlFG4_Sc6xY)

I used a pure numpy functions which helped me to calculate the IK fast, tf.transformations is based on numpy, the functions message_from_transform, ik_calculator.FK_calculator, ik_calculator.joint_states_callback was for ploting error purpose until now I didn't implement it .

I used inverse of the correction matrix to change the input from urdf to DH representation .

the calculations of angles are as show in Kinematic Analysis section .




