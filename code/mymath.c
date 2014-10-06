#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mymath.h"

double_triple add_double_triple(double_triple a, double_triple b){
	double_triple result;
	result.x=(a).x+(b).x;
	result.y=(a).y+(b).y;
	result.z=(a).z+(b).z;
	return result;
}

double_triple scalar_multiply_double_triple(double_triple a, double s){
	double_triple result;
	result.x=s*a.x;
	result.y=s*a.y;
	result.z=s*a.z;
	return result;
}

double dot_product(double_triple a, double_triple b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

double norm(double_triple a){
	return sqrt(dot_product(a, a));
}

double_triple rand_unit_cube(){
	double_triple result;
	result.x=2*rand_double-1;
	result.y=2*rand_double-1;
	result.z=2*rand_double-1;
	return result;
}

double_triple rand_unit_ball(){
	double_triple result;
	double r2=1;
	while(r2>=1){
		result.x=2*rand_double-1;
		result.y=2*rand_double-1;
		result.z=2*rand_double-1;
		r2=dot_product(result, result);
	}
	return result;
}

double_triple rand_unit_sphere(){
	double_triple result;
	double r2=1;
	while(r2>=1){
		result.x=2*rand_double-1;
		result.y=2*rand_double-1;
		result.z=2*rand_double-1;
		r2=dot_product(result, result);
	}
	r2=1/(sqrt(r2));
	result.x*=r2;
	result.y*=r2;
	result.z*=r2;
	return result;
}

double_triple subtract_double_triple(double_triple a, double_triple b){
	double_triple result;
	result.x=a.x-b.x;
	result.y=a.y-b.y;
	result.z=a.z-b.z;
	return result;
}

void recenter_double_triple(double_triple *pa, double_triple b){
	recenter(((*pa).x), ((b).x));
	recenter(((*pa).y), ((b).y));
	recenter(((*pa).z), ((b).z));
}

void fmod_double_triple(double_triple *pa, double_triple b){
	fmod(((*pa).x), ((b).x));
	fmod(((*pa).y), ((b).y));
	fmod(((*pa).z), ((b).z));
}

void normalize(double_triple *pvec){
	double mynorm=norm(*pvec);
	(*pvec).x/=mynorm;
	(*pvec).y/=mynorm;
	(*pvec).z/=mynorm;
}

double_triple normed(double_triple vec){
	double_triple result=(vec);
    double oneovernorm=1./norm(vec);
    result.x*=oneovernorm;
    result.y*=oneovernorm;
    result.z*=oneovernorm;
    return result;
}

double_triple cross_product(double_triple a, double_triple b){
	double_triple result;
	result.x=a.y*b.z-a.z*b.y;
	result.y=a.z*b.x-a.x*b.z;
	result.z=a.x*b.y-a.y*b.x;
	return result;
}

double largest_cubic_root(double a3, double a2, double a1, double a0){
    
    // from http://www-old.me.gatech.edu/energy/andy_phd/appA.htm
    // except he was missing a 1/3 in definition of phi
    // corrected by looking at http://home.pipeline.com/~hbaker1/cubic3realroots.htm
    
    double p=a2/a3;
    double q=a1/a3;
    double r=a0/a3;
    double A=(3*q-p*p)/3;
    double B=(2*p*p*p-9*p*q+27*r)/27;
    double D=A*A*A/27+B*B/4;
    double M, N, y1, y2, y3, phi, x1, x2, x3, biggest, M3, N3, arg;
    //printf("finding root of %f, %f, %f, %f\n", a3, a2, a1, a0);
    //printf("D=%f\n", D);
    if(D>0){                        //  one real root
        M3=-B/2+sqrt(D);
        if(M3<0){
            M=-pow(-M3, 1./3.);
            //printf("M=%f^(1/3)=%f\n", M3, M);
        }
        else{
            M=pow(M3, 1./3.);
            //printf("M=%f^(1/3)=%f\n", M3, M);
        }
        N3=-B/2-sqrt(D);
        if(N3<0){
            N=-pow(-N3, 1./3.);
            //printf("N=%f^(1/3)=%f\n", N3, N);
        }
        else{
            N=pow(N3, 1./3.);
            //printf("N=%f^(1/3)=%f\n", N3, N);
        }
        y1=M+N;
        x1=y1-p/3;
        biggest=x1;
    }
    else{                           //  three real roots
        if(B>0){
            arg=B*B/4/(-A*A*A/27);
            //printf("arg=%f\n", arg);
            phi=acos(-sqrt(arg))/3;
        }
        else{
            arg=B*B/4/(-A*A*A/27);
            //printf("arg=%f\n", arg);
            phi=acos(sqrt(arg))/3;
        }
        //printf("A=%f\n", A);
        y1=2*sqrt(-A/3)*cos(phi);
        y2=2*sqrt(-A/3)*cos(phi+2./3.*M_PI);
        y3=2*sqrt(-A/3)*cos(phi+4./3.*M_PI);
        x1=y1-p/3;
        x2=y2-p/3;
        x3=y3-p/3;
        //printf("roots=%f, %f, %f\n", x1, x2, x3);
        if(x1>x2){
            if(x1>x3) biggest=x1;
            else biggest=x3;
        }
        else{
            if(x2>x3) biggest=x2;
            else biggest=x3;
        }
    }
    //printf("root %f\n", biggest);
    return biggest;
}

double quartic_global_minimum(double p4, double p3, double p2, double p1){
    double a3=4*p4;
    double a2=3*p3;
    double a1=2*p2;
    double a0=p1;

    // solve cubic
    
    // from http://www-old.me.gatech.edu/energy/andy_phd/appA.htm
    // except he was missing a 1/3 in definition of phi
    // corrected by looking at http://home.pipeline.com/~hbaker1/cubic3realroots.htm
    
    double p=a2/a3;
    double q=a1/a3;
    double r=a0/a3;
    double A=(3*q-p*p)/3;
    double B=(2*p*p*p-9*p*q+27*r)/27;
    double D=A*A*A/27+B*B/4;
    double M, N, y1, y2, y3, phi, x1, x2, x3, lowest, M3, N3, arg;
    //printf("finding root of %f, %f, %f, %f\n", a3, a2, a1, a0);
    //printf("D=%f\n", D);
    if(D>0){                        //  one real root
        M3=-B/2+sqrt(D);
        if(M3<0){
            M=-pow(-M3, 1./3.);
            //printf("M=%f^(1/3)=%f\n", M3, M);
        }
        else{
            M=pow(M3, 1./3.);
            //printf("M=%f^(1/3)=%f\n", M3, M);
        }
        N3=-B/2-sqrt(D);
        if(N3<0){
            N=-pow(-N3, 1./3.);
            //printf("N=%f^(1/3)=%f\n", N3, N);
        }
        else{
            N=pow(N3, 1./3.);
            //printf("N=%f^(1/3)=%f\n", N3, N);
        }
        y1=M+N;
        x1=y1-p/3;
        lowest=x1;
    }
    else{                           //  three real roots
        if(B>0){
            arg=B*B/4/(-A*A*A/27);
            //printf("arg=%f\n", arg);
            phi=acos(-sqrt(arg))/3;
        }
        else{
            arg=B*B/4/(-A*A*A/27);
            //printf("arg=%f\n", arg);
            phi=acos(sqrt(arg))/3;
        }
        //printf("A=%f\n", A);
        y1=2*sqrt(-A/3)*cos(phi);
        y2=2*sqrt(-A/3)*cos(phi+2./3.*M_PI);
        y3=2*sqrt(-A/3)*cos(phi+4./3.*M_PI);
        x1=y1-p/3;
        x2=y2-p/3;
        x3=y3-p/3;
        //printf("roots=%f, %f, %f\n", x1, x2, x3);
        if(p4*pow(x1, 4)+p3*pow(x1, 3)+p2*pow(x1, 2)+p1*x1<p4*pow(x2, 4)+p3*pow(x2, 3)+p2*pow(x2, 2)+p1*x2){
            if(p4*pow(x1, 4)+p3*pow(x1, 3)+p2*pow(x1, 2)+p1*x1<p4*pow(x3, 4)+p3*pow(x3, 3)+p2*pow(x3, 2)+p1*x3){
                lowest=x1;
            }
            else lowest=x3;
        }
        else{
            if(p4*pow(x2, 4)+p3*pow(x2, 3)+p2*pow(x2, 2)+p1*x2<p4*pow(x3, 4)+p3*pow(x3, 3)+p2*pow(x3, 2)+p1*x3){
                lowest=x2;
            }
            else lowest=x3;
        }
    }
    //printf("root %f\n", biggest);
    return lowest;
}

void matrix_multiply(double **a, double **b, double **c){
    c[0][0]=a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0];
    c[0][1]=a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1];
    c[0][2]=a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2];
    c[1][0]=a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0];
    c[1][1]=a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1];
    c[1][2]=a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2];
    c[2][0]=a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0];
    c[2][1]=a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1];
    c[2][2]=a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2];
}

void forward_matrix(double angle, double_triple axis_vector, double **forward){
	double costheta=cos(angle);
	double oneminuscos=1.0-cos(angle);
	double sintheta=sin(angle);
	forward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	forward[0][1]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	forward[0][2]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	forward[1][0]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	forward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	forward[1][2]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	forward[2][0]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	forward[2][1]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	forward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
}

void forward_and_backward_matrix(double angle, double_triple axis_vector, double **forward, double **backward){
	double costheta=cos(angle);
	double oneminuscos=1.0-cos(angle);
	double sintheta=sin(angle);
	forward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	forward[0][1]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	forward[0][2]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	forward[1][0]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	forward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	forward[1][2]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	forward[2][0]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	forward[2][1]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	forward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
	backward[0][0]=costheta+axis_vector.x*axis_vector.x*oneminuscos;
	backward[0][1]=axis_vector.x*axis_vector.y*oneminuscos+axis_vector.z*sintheta;
	backward[0][2]=axis_vector.x*axis_vector.z*oneminuscos-axis_vector.y*sintheta;
	backward[1][0]=axis_vector.x*axis_vector.y*oneminuscos-axis_vector.z*sintheta;
	backward[1][1]=costheta+axis_vector.y*axis_vector.y*oneminuscos;
	backward[1][2]=axis_vector.y*axis_vector.z*oneminuscos+axis_vector.x*sintheta;
	backward[2][0]=axis_vector.x*axis_vector.z*oneminuscos+axis_vector.y*sintheta;
	backward[2][1]=axis_vector.y*axis_vector.z*oneminuscos-axis_vector.x*sintheta;
	backward[2][2]=costheta+axis_vector.z*axis_vector.z*oneminuscos;
}

double_triple rotate_by_matrix(double_triple start, double **matrix){
	double_triple result;
	result.x=start.x*matrix[0][0]+start.y*matrix[0][1]+start.z*matrix[0][2];
	result.y=start.x*matrix[1][0]+start.y*matrix[1][1]+start.z*matrix[1][2];
	result.z=start.x*matrix[2][0]+start.y*matrix[2][1]+start.z*matrix[2][2];
	return result;
}

void rand_unit_quaternion(double_quadruple *punit){
	double mag=2.0;
	double_quadruple temp;
	while(mag>1.0){
		temp.q0=1.0-2.0*rand_double;
		temp.q1=1.0-2.0*rand_double;
		temp.q2=1.0-2.0*rand_double;
		temp.q3=1.0-3.0*rand_double;
		mag=temp.q0*temp.q0+temp.q1*temp.q1+temp.q2*temp.q2+temp.q3*temp.q3;
	}
	mag=1.0/sqrt(mag);
	(*punit).q0=temp.q0*mag;
	(*punit).q1=temp.q1*mag;
	(*punit).q2=temp.q2*mag;
	(*punit).q3=temp.q3*mag;
}

void normalize_quaternion(double_quadruple *pquaternion){
	double factor=((*pquaternion).q0)*((*pquaternion).q0)+((*pquaternion).q1)*((*pquaternion).q1)+((*pquaternion).q2)*((*pquaternion).q2)+((*pquaternion).q3)*((*pquaternion).q3);
	factor=1.0/sqrt(factor);
	((*pquaternion).q0)*=factor;
	((*pquaternion).q1)*=factor;
	((*pquaternion).q2)*=factor;
	((*pquaternion).q3)*=factor;
}

double_quadruple quaternion_product(double_quadruple a, double_quadruple b){
	double_quadruple result;
	result.q0=a.q0*b.q0-a.q1*b.q1-a.q2*b.q2-a.q3*b.q3;
	result.q1=a.q0*b.q1+a.q1*b.q0+a.q2*b.q3-a.q3*b.q2;
	result.q2=a.q0*b.q2-a.q1*b.q3+a.q2*b.q0+a.q3*b.q1;
	result.q3=a.q0*b.q3+a.q1*b.q2-a.q2*b.q1+a.q3*b.q0;
	return result;
}

double_quadruple quaternion_reverse_product(double_quadruple a, double_quadruple b){
	double_quadruple result;
	b.q1=-b.q1;
	b.q2=-b.q2;
	b.q3=-b.q3;
	result.q0=a.q0*b.q0-a.q1*b.q1-a.q2*b.q2-a.q3*b.q3;
	result.q1=a.q0*b.q1+a.q1*b.q0+a.q2*b.q3-a.q3*b.q2;
	result.q2=a.q0*b.q2-a.q1*b.q3+a.q2*b.q0+a.q3*b.q1;
	result.q3=a.q0*b.q3+a.q1*b.q2-a.q2*b.q1+a.q3*b.q0;
	return result;
}

double_triple rotate_by_quaternion(double_triple vector, double_quadruple a){
	double_triple result;
	result.x=(1.0-2.0*(a.q2*a.q2+a.q3*a.q3))*vector.x+2.0*(a.q1*a.q2+a.q0*a.q3)*vector.y+2.0*(a.q1*a.q3-a.q0*a.q2)*vector.z;
	result.y=2.0*(a.q1*a.q2-a.q0*a.q3)*vector.x+(1.0-2.0*(a.q1*a.q1+a.q3*a.q3))*vector.y+2.0*(a.q2*a.q3+a.q0*a.q1)*vector.z;
	result.z=2.0*(a.q1*a.q3+a.q0*a.q2)*vector.x+2.0*(a.q2*a.q3-a.q0*a.q1)*vector.y+(1.0-2.0*(a.q1*a.q1+a.q2*a.q2))*vector.z;
	return result;
}

double quaternion_dot_product(double_quadruple a, double_quadruple b){
	double_triple avec, bvec;
	avec.x=(1.0-2.0*(a.q2*a.q2+a.q3*a.q3));
	avec.y=2.0*(a.q1*a.q2-a.q0*a.q3);
	avec.z=2.0*(a.q1*a.q3+a.q0*a.q2);
	bvec.x=(1.0-2.0*(b.q2*b.q2+b.q3*b.q3));
	bvec.y=2.0*(b.q1*b.q2-b.q0*b.q3);
	bvec.z=2.0*(b.q1*b.q3+b.q0*b.q2);
	return dot_product(avec, bvec);
}

double_triple rotate_by_quaternion_reverse(double_triple vector, double_quadruple a){
    
    //	rotate by inverse of a: (a.q0, a.q1, a.q2, a.q3) -> (a.q0, -a.q1, -a.q2, -a.q3)
	
    double_triple result;
	result.x=(1.0-2.0*(a.q2*a.q2+a.q3*a.q3))*vector.x+2.0*(a.q1*a.q2-a.q0*a.q3)*vector.y+2.0*(a.q1*a.q3+a.q0*a.q2)*vector.z;
	result.y=2.0*(a.q1*a.q2+a.q0*a.q3)*vector.x+(1.0-2.0*(a.q1*a.q1+a.q3*a.q3))*vector.y+2.0*(a.q2*a.q3-a.q0*a.q1)*vector.z;
	result.z=2.0*(a.q1*a.q3+a.q0*a.q2)*vector.x+2.0*(a.q2*a.q3+a.q0*a.q1)*vector.y+(1.0-2.0*(a.q1*a.q1+a.q2*a.q2))*vector.z;
	return result;
}

void invert_3x3_matrix(double **a, double **b){
    int i, j;
    double det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){            
            b[i][j]=(a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3])/det;
        }
    }
}

