#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"
gsl_rng *tau;

double norma(double x,double y)
{
double norma;
norma=sqrt((x*x)+(y*y));
return norma;                                       
}

double dist(double xi, double yi, double xj, double yj)
{
double a,b,d,e;
a=xi-xj;
b=yi-yj;		
if(fabs(a)>2.0)
	{
	if(xj>2.0)
	 d=xj-4; 
	 else
	 d=xj+4;
	a=xi-d;
	}
if(fabs(b)>2.0)
	{
	if(yj>2.0)
	 e=yj-4; 
	 else
	 e=yj+4;
	b=yi-e;
	}		
return norma(a,b);
}


int main()
{
double x[16],y[16],vx[16],vy[16],ax[16],ay[16],wx[16],wy[16],ptransf[16],h=0.0001,n,m,E,K,U,r[16][16],a,b,c=2.0,e,d,temp[1000];
double temperatura,P,Ptransf=0.0,presion,fluct,x0,y0,scm,Tc[1000],TempCrit;
int i,j,k,l,q,t=0,control=0,g=0;
long int seed=447263571;
extern gsl_rng *tau;
FILE *f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9,*f10;
tau=gsl_rng_alloc(gsl_rng_taus);
gsl_rng_set(tau,seed);
f1=fopen("posiciones.txt","w");
f2=fopen("comprobar.txt","w");
f3=fopen("velocidad0.txt","w");
f4=fopen("energia.txt","w");
f5=fopen("a_ini.txt","w");
f6=fopen("r_ini.txt","w");
f7=fopen("temperatura.txt","w");
f8=fopen("Fluctuacionr.txt","w");
f9=fopen("SepCuadMed.txt","w");
f10=fopen("TempCrit.txt","w");
                                           //Posiciones iniciales               
                       //Ordenados
if(control==0)
{
x[0]=1.25*4.0/10.0;
y[0]=8.75*4.0/10.0;
x[1]=3.75*4.0/10.0;
y[1]=8.75*4.0/10.0;
x[2]=6.25*4.0/10.0;
y[2]=8.75*4.0/10.0;
x[3]=8.75*4.0/10.0;
y[3]=8.75*4.0/10.0;
x[4]=1.25*4.0/10.0;
y[4]=6.25*4.0/10.0;
x[5]=3.75*4.0/10.0;
y[5]=6.25*4.0/10.0;
x[6]=6.25*4.0/10.0;
x0=x[6];
y[6]=6.25*4.0/10.0;
y0=y[6];
x[7]=8.75*4.0/10.0;
y[7]=6.25*4.0/10.0;
x[8]=1.25*4.0/10.0;
y[8]=3.75*4.0/10.0;
x[9]=3.75*4.0/10.0;
y[9]=3.75*4.0/10.0;
x[10]=6.25*4.0/10.0;
y[10]=3.75*4.0/10.0;
x[11]=8.75*4.0/10.0;
y[11]=3.75*4.0/10.0;
x[12]=1.25*4.0/10.0;
y[12]=1.25*4.0/10.0;
x[13]=3.75*4.0/10.0;
y[13]=1.25*4.0/10.0;
x[14]=6.25*4.0/10.0;
y[14]=1.25*4.0/10.0;
x[15]=8.75*4.0/10.0;
y[15]=1.25*4.0/10.0;
}
                    //Aleatorial
if (control==1)
{
for(i=0;i<16;i++)
{
n=gsl_rng_uniform(tau);
m=gsl_rng_uniform(tau);
x[i]=4.0*n;
y[i]=4.0*m;
}
}




                                           //Velocidades iniciales
for(i=0;i<16;i++)
	{
	n=gsl_rng_uniform(tau);
	vx[i]=0.3*cos(n*2.0*3.141592);
	vy[i]=0.3*sin(n*2.0*3.141592);
	fprintf(f3,"%lf\t%lf\t%lf\n",vx[i],vy[i],sqrt(vx[i]*vx[i]+vy[i]*vy[i]));
	}
	               //Setear a[] y w[] a 0
for(j=0;j<16;j++)
	{
	ax[j]=0.0;
	ay[j]=0.0;
	wx[j]=0.0;
	wy[j]=0.0;
	}
		     //Calcular a[] en t a partir de (x,y) (con cond.period.)		   
for(i=0;i<16;i++)
	for(j=0;j<16;j++)
		{
		if(i!=j)
		{
a=x[i]-x[j];
b=y[i]-y[j];		
if(fabs(a)>c)
	{
	if(x[j]>c)
	 d=x[j]-4; 
	 else
	 d=x[j]+4;
	a=x[i]-d;
	}
if(fabs(b)>c)
	{
	if(y[j]>c)
	 e=y[j]-4; 
	 else
	 e=y[j]+4;
	b=y[i]-e;
	}		
r[i][j]=norma(a,b);
ax[i]=ax[i]-24*a*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));
ay[i]=ay[i]-24*b*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));		
		}	
		}
                                  //Plot aceleraciones iniciales
for(j=0;j<16;j++)
	{
	fprintf(f5,"%lf\t%lf\n",ax[j],ay[j]);
	}		
		
for(q=0;q<800000;q++)
{	
for (i=0;i<16;i++)
	ptransf[i]=0;
		
		    //Calcular posiciones en t+h
for(j=0;j<16;j++)
	{
	x[j]=x[j]+(h*vx[j])+(0.5*h*h*ax[j]);
	y[j]=y[j]+(h*vy[j])+(0.5*h*h*ay[j]);
	if(x[j]>4)
	{
	x[j]=(x[j]-4);
	ptransf[j]=2*vx[j];
	}
	if(x[j]<0)
	x[j]=x[j]+4;
	if(y[j]>4)
	y[j]=y[j]-4;
	if(y[j]<0)
	y[j]=y[j]+4;
	}	
                       		  
	
	               //Plotear posiciones
if(q%200==0)
{	               
for(j=0;j<15;j++)
	{
	fprintf(f1,"%lf\t%lf\t",x[j],y[j]);
	}
fprintf(f1,"%lf\t%lf\t\n",x[15],y[15]);
}
	  
	             //Calcular las w[] en t
for(i=0;i<16;i++)
	{
	wx[i]=vx[i]+(0.5*h*ax[i]);
	wy[i]=vy[i]+(0.5*h*ay[i]);
	}	             
	           //Calcular aceleraciones en t+h a partir de (x,y) y en t+h
for(l=0;l<16;l++)
	{
	ax[l]=0.0;
	ay[l]=0.0;
	}	
	           
for(i=0;i<16;i++)
	for(j=0;j<16;j++)
		{
		if(i!=j)
		{
a=x[i]-x[j];
b=y[i]-y[j];		
if(fabs(a)>c)
	{
	if(x[j]>c)
	 d=x[j]-4; 
	 else
	 d=x[j]+4;
	a=x[i]-d;
	}
if(fabs(b)>c)
	{
	if(y[j]>c)
	 e=y[j]-4; 
	 else
	 e=y[j]+4;
	b=y[i]-e;
	}		
r[i][j]=norma(a,b);
ax[i]=ax[i]-24*a*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));
ay[i]=ay[i]-24*b*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));		            
		}	
		}
		        //Calc v=w+(h²/2)a
for(i=0;i<16;i++)
	{
	vx[i]=wx[i]+(0.5*h*h*ax[i]);
	vy[i]=wy[i]+(0.5*h*h*ay[i]);
	if(q%10000==0)
		{
		vx[i]=1.1*vx[i];
		vy[i]=1.1*vy[i];
		}
	}
	
	                    //Calculo E
if(q%50==0)
{
E=0.0;
K=0.0;
U=0.0;
	          
for(i=0;i<16;i++)
	{
	K=K+(0.5*pow(norma(vx[i],vy[i]),2));
	}
for(j=0;j<16;j++)
	{
	for(k=0;k<16;k++)
	{
	if(k!=j)
		U=U+4*((pow(r[j][k],-12))-(pow(r[j][k],-6)));
	}
	}	
	
E=K+U;	
	              //Escribir en fichero E
fprintf(f4,"%i\t%lf\t%lf\t%lf\n",q,E,K,U);				
}

			                
for(j=0;j<16;j++)
 {
 if((x[j]>(50))||(y[j]>(50))||(x[j]<(-50))||(y[j]<(-50)))
	fprintf(f2,"error part %i en el paso %i\n",j,q);	
}


                       //Calculo de fluctuaciones en posicion
fluct=0.0;
fluct=pow((dist(x[6],y[6],x0,y0)),2);                       
if(q%100==0)
fprintf(f8,"%i\t%lf\n",q,fluct);                       
                                //Calculo de separacion cuadratica media
scm=0.0;
scm=pow((dist(x[6],y[6],x[7],y[7])),2);
if(q%100==0)
fprintf(f9,"%i\t%lf\n",q,scm); 

                             //Calculo de temp Crítica
if((q>186500)&&(q<187500))
{
Tc[g]=0.0;
for(i=0;i<16;i++)
{
temp[g]=temp[g]+(0.5*((vx[i]*vx[i])+(vy[i]*vy[i])));
}
fprintf(f10,"%lf\n",temp[g]);
g++;
}

}                                   //187000          240000



fclose(f10);
fclose(f9);                                          
fclose(f8);
fclose(f7);	
fclose(f6);
fclose(f5);
fclose(f4);		
fclose(f3);
fclose(f2);
fclose(f1);
//printf("Temperatura %lf\n",temperatura);
//printf("Presion %lf\n",presion);
return 0;
}
          
