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


int main()
{
double x[20],y[20],vx[20],vy[20],ax[20],ay[20],wx[20],wy[20],ptransf[20],h=0.002,n,m,g,E,K,U,r[20][20],a,b,c=5.0,e,d,temp[1000];
double temperatura,P,Ptransf=0.0,presion;
int i,j,k,l,q,t=0,cont=0;
long int seed=447263571;
extern gsl_rng *tau;
FILE *f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9;
tau=gsl_rng_alloc(gsl_rng_taus);
gsl_rng_set(tau,seed);
f1=fopen("posiciones.txt","w");
f2=fopen("comprobar.txt","w");
f3=fopen("velocidad0.txt","w");
f4=fopen("energia.txt","w");
f5=fopen("a_ini.txt","w");
f6=fopen("r_ini.txt","w");
f7=fopen("histograma.txt","w");
f8=fopen("presion.txt","w");
f9=fopen("temperatura.txt","w");
                                           //Posiciones iniciales                                   
                       
                       //Aleatorio.
for(i=0;i<20;i++)
	{
	n=gsl_rng_uniform(tau);
	m=gsl_rng_uniform(tau);
	x[i]=10.0*n;
	y[i]=10.0*m;
	fprintf(f6,"%lf\t%lf\n",x[i],y[i]);
	}


                                           //Velocidades iniciales
for(i=0;i<20;i++)
	{
	g=gsl_rng_uniform(tau);
	vx[i]=12.0*cos(n*2.0*3.141592);
	vy[i]=12.0*sin(n*2.0*3.141592);
	fprintf(f3,"%lf\t%lf\t%lf\n",vx[i],vy[i],sqrt(vx[i]*vx[i]+vy[i]*vy[i]));
	}
	               //Setear a[] y w[] a 0
for(j=0;j<20;j++)
	{
	ax[j]=0.0;
	ay[j]=0.0;
	wx[j]=0.0;
	wy[j]=0.0;
	}
		     //Calcular a[] en t a partir de (x,y) (con cond.period.)		   
for(i=0;i<20;i++)
	for(j=0;j<20;j++)
		{
		if(i!=j)
		{
a=x[i]-x[j];
b=y[i]-y[j];		
if(fabs(a)>c)
	{
	if(x[j]>c)
	 d=x[j]-10; 
	 else
	 d=x[j]+10;
	a=x[i]-d;
	}
if(fabs(b)>c)
	{
	if(y[j]>c)
	 e=y[j]-10; 
	 else
	 e=y[j]+10;
	b=y[i]-e;
	}		
r[i][j]=norma(a,b);
ax[i]=ax[i]-24*a*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));
ay[i]=ay[i]-24*b*((-2*pow(r[i][j],-14))+(pow(r[i][j],-8)));		
		}	
		}
                                  //Plot aceleraciones iniciales
for(j=0;j<20;j++)
	{
	fprintf(f5,"%lf\t%lf\n",ax[j],ay[j]);
	}		
		
for(q=0;q<160000;q++)
{	
for (i=0;i<20;i++)
	ptransf[i]=0;
		
		    //Calcular posiciones en t+h
for(j=0;j<20;j++)
	{
	x[j]=x[j]+(h*vx[j])+(0.5*h*h*ax[j]);
	y[j]=y[j]+(h*vy[j])+(0.5*h*h*ay[j]);
	if(x[j]>10)
	{
	x[j]=(x[j]-10);
	ptransf[j]=2*vx[j];                    //momento transferido por atomo j en esta iteracion
	}
	if(x[j]<0)
	x[j]=x[j]+10;
	if(y[j]>10)
	y[j]=y[j]-10;
	if(y[j]<0)
	y[j]=y[j]+10;
	}	
                       		  
	
	               //Plotear posiciones
if(q%100==0)
{	               
for(j=0;j<19;j++)
	{
	fprintf(f1,"%lf\t%lf\t",x[j],y[j]);
	}
fprintf(f1,"%lf\t%lf\t\n",x[19],y[19]);
}
	  
	             //Calcular las w[] en t
for(i=0;i<20;i++)
	{
	wx[i]=vx[i]+(0.5*h*ax[i]);
	wy[i]=vy[i]+(0.5*h*ay[i]);
	}	             
	           //Calcular aceleraciones en t+h a partir de (x,y) y en t+h
for(l=0;l<20;l++)
	{
	ax[l]=0.0;
	ay[l]=0.0;
	}	
	           
for(i=0;i<20;i++)
	for(j=0;j<20;j++)
		{
		if(i!=j)
		{
a=x[i]-x[j];
b=y[i]-y[j];		
if(fabs(a)>c)
	{
	if(x[j]>c)
	 d=x[j]-10; 
	 else
	 d=x[j]+10;
	a=x[i]-d;
	}
if(fabs(b)>c)
	{
	if(y[j]>c)
	 e=y[j]-10; 
	 else
	 e=y[j]+10;
	b=y[i]-e;
	}		
r[i][j]=norma(a,b);
ax[i]=ax[i]-24.0*a*((-2.0*pow(r[i][j],-14))+(pow(r[i][j],-8)));
ay[i]=ay[i]-24.0*b*((-2.0*pow(r[i][j],-14))+(pow(r[i][j],-8)));		            
		}	
		}
		        //Calc v=w+(h²/2)a
for(i=0;i<20;i++)
	{
	vx[i]=wx[i]+(0.5*h*h*ax[i]);
	vy[i]=wy[i]+(0.5*h*h*ay[i]);
	}
	                    //Calculo E
if(q%10==0)
{
E=0.0;
K=0.0;
U=0.0;
	          
for(i=1;i<20;i++)
	{
	K=K+(0.5*pow(norma(vx[i],vy[i]),2));
	}
for(j=1;j<16;j++)
	{
	for(k=1;k<16;k++)
	{
	if(k!=j)
		U=U+4.0*((pow(r[j][k],-12))-(pow(r[j][k],-6)));
	}
	}		
	
E=K+U;	
	              //Escribir en fichero E
fprintf(f4,"%i\t%lf\t%lf\t%lf\n",q,E,K,U);				
}

	            //Comprobar si alguna se ha ido a + de 50          
for(j=0;j<20;j++)
 {
 if((x[j]>(50))||(y[j]>(50))||(x[j]<(-50))||(y[j]<(-50)))
	fprintf(f2,"error part %i en el paso %i\n",j,q);	
}

                         //Calculo de Temperatura
if((q>10000)&&(q<11000))
{
temp[t]=0;
for(i=0;i<20;i++)
{
temp[t]=temp[t]+((0.5/20.0)*((vx[i]*vx[i])+(vy[i]*vy[i])));

}
t++;

if(q==10750)         //Hist. vielocidades
	{
	for(i=0;i<20;i++)
	fprintf(f7,"%lf\n",norma(vx[i],vy[i]));
	}

}
temperatura=0.0;
for(j=0;j<1000;j++)
{
temperatura=temperatura+(0.001*temp[j]);          
if((q%100==0)&&(q>50000))
fprintf(f9,"%lf\n",temp[j]);                    //Plotear temperat
}
                                 //Calculo del momento total en esta iteración
for(i=0;i<20;i++)
Ptransf=Ptransf+ptransf[i];


P=(Ptransf/(q*h*10));	       //Calculo de la presión en esta iteración

if((q%50==0)&&(q!=0))	         //Plotear presiones pa sacar desvest.
fprintf(f8,"%lf\n",P);	  
            
if((q<10000)&&(q>8000))            //promediar presión entre 8000 y 10000
presion=presion+(P/2000);



}
fclose(f9);
fclose(f8);
fclose(f7);	
fclose(f6);
fclose(f5);
fclose(f4);		
fclose(f3);
fclose(f2);
fclose(f1);
printf("Temperatura %lf\n",temperatura);
printf("Presion %lf\n",presion);
return 0;
}
