//==============================================================================
//
// Title:       Genetic Algorithm routines
// Purpose:     Algorithms for performing genetic optimization
//
// Created on:  14-09-2012 by Bas van der Geer
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================


#include "toolbox.h"
#include <userint.h>
#include <ansi_c.h>
#include "SLM_genetic.h" 

static unsigned int gPopSize;
static unsigned int gVarSize;
static double scale=0.6 ;
static double rho=0.5 ;
static unsigned int toursize=3 ; 
static dFpd gQualityF;


struct species
{
	tGeneticVar *x ;
	double f ;
} *gPopulation ;


struct species *get_good(void)
{
	int ibest ;
	double fbest=1.0e100 ;

	for(unsigned int i=0 ; i<toursize ; i++)
	{
		unsigned int j = rand()%gPopSize ;
		if( gPopulation[j].f < fbest )
		{
			fbest = gPopulation[j].f ;
			ibest=j ;
		}
	}
	//printf("  good: %u, %.5g \n", ibest, fbest);
	return &gPopulation[ibest] ;
}

struct species *get_bad(void)
{
	int iworst ;
	double fbest=0 ;

	for(unsigned int i=0 ; i<toursize ; i++)
	{
		unsigned int j = rand()%gPopSize ;
		if( gPopulation[j].f > fbest )
		{
			fbest = gPopulation[j].f ;
			iworst=j ;
		}
	}
	//printf("  bad:  %u, %.5g \n", iworst, fbest);
	return &gPopulation[iworst] ;
}


/*
double f(tGeneticVar * restrict phase)
{
	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI);
	
	// compose the field in the SLM plane from the amplitude and the supplied phase
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
	{
		gFFTin[l * gXsize + k][0] = gInputAmplitude[l * gXsize + k] * cos(((double) phase[l * gXsize + k]) / convf);
		gFFTin[l * gXsize + k][1] = gInputAmplitude[l * gXsize + k] * sin(((double) phase[l * gXsize + k]) / convf);
	}

	// propagate the light field forwards to obtain the (simulated) intensity in the focal plane
	fftw_execute(gFFTplan_Gg);	
	
	// measure the total intensity in the signal window
	double winint = 0.0;
	for (int k = gWindowX; k < gWindowX + gWindowXsize; k++)
	for (int l = gWindowX; l < gWindowX + gWindowXsize; l++)
		winint += ( gFFTout[l * gXsize + k][0] * gFFTout[l * gXsize + k][0] 
			      + gFFTout[l * gXsize + k][1] * gFFTout[l * gXsize + k][1] );
	
	double sigconv = gSigInt / winint;
	
	
	// compute the quality measure (lower is better)
	double quality = 0.0;
	
	// sum only over the pixels in the specified window area
	for (int k = gWindowX; k < gWindowX + gWindowXsize; k++)
	for (int l = gWindowY; l < gWindowY + gWindowYsize; l++)
	{
		// compute the deviation in amplitude
		double abssig = sqrt( gFFTout[l * gXsize + k][0] * gFFTout[l * gXsize + k][0] 
			                + gFFTout[l * gXsize + k][1] * gFFTout[l * gXsize + k][1] ) * sqrt(sigconv);
		double deviation = fabs(gSignal[l * gXsize + k] - abssig);
		
		// add the square of the deviation to the 'quality'
		quality += deviation * deviation;
	}
	
	return quality / (gXsize * gYsize);
}*/


tGeneticVar combine(tGeneticVar u, tGeneticVar v, tGeneticVar w, double scale)
{
	double t = w + scale * (u - v) ;

	//if(t<0) t=0 ;
	//if(t>255) t=255 ;

	return (tGeneticVar) t ;
}


tGeneticVar crossover(tGeneticVar i, tGeneticVar s,double rho)
{
	if( Random(0.0, 1.0) < rho )
		return i ;
	else
		return s ;
}


void GeneticInitialize(int popsize, unsigned int numgenes, tGeneticVar maxvariance, dFpd f)
{
	// was everything already initialised?
	if (gPopulation != NULL)
	{
		// free all individual pointers
		for(unsigned int i = 0; i < gPopSize + 1; i++)
			free(gPopulation[i].x);
		
		// free the entire population pointer
		free(gPopulation);
	}
	
	// copy (new) population size and number of genes (gVarSize)
	gPopSize = popsize;
	gVarSize = numgenes;
	
	// set the new quality measure function pointer
	gQualityF = f;

	// Initialize new population
	// (one extra element, used as a temporary)
	gPopulation = (struct species *)malloc((gPopSize + 1) * sizeof(struct species)) ;

	// fill the population
	for(unsigned int i = 0; i < gPopSize + 1; i++)
	{
		gPopulation[i].x = (tGeneticVar *) malloc(gVarSize * sizeof(tGeneticVar)) ;
		for(unsigned int j = 0; j < gVarSize; j++) 
			gPopulation[i].x[j] = Random(-1.0, 1.0) * maxvariance;
		gPopulation[i].f = f(gPopulation[i].x) ;
	}
}



int GeneticOptimize(int NumEval, int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6, int debuggraph1)
{
	// Our best attempt so far
	unsigned int ibest ;
	double fbest=DBL_MAX ;
	
	// Get best
	for(unsigned int i=0 ; i<gPopSize ; i++)
		if( gPopulation[i].f < fbest )
	{
		fbest=gPopulation[i].f ;
		ibest=i ;
	}

	// Run DE
	unsigned int count = 0;
	double* fbestit = (double*) calloc(NumEval, sizeof(double));
	double* fnew    = (double*) calloc(NumEval, sizeof(double));
	double* fworstit    = (double*) calloc(NumEval, sizeof(double));
	double* favg    = (double*) calloc(NumEval, sizeof(double));
	while(count < NumEval)
	{
		// store the current best score
		fbestit[count] = fbest;
		
		// Counters
		count++ ;

		// Four random parents
		struct species *xi = get_bad();  // &gPopulation[rand()%gPopSize] ; // substitution candidate
		struct species *xu = get_good(); // &gPopulation[rand()%gPopSize] ;
		struct species *xv = get_good(); // &gPopulation[rand()%gPopSize] ;
		struct species *xw = get_good(); // &gPopulation[rand()%gPopSize] ;

		struct species *xs = &gPopulation[gPopSize] ; // Our temp variable

		// Combine
		for(unsigned int j=0; j < gVarSize ; j++)
			xs->x[j] = combine(xu->x[j],xv->x[j],xw->x[j],scale) ;

		// Crossover
		for(unsigned int j=0 ; j<gVarSize ; j++)
			xs->x[j] = crossover(xi->x[j],xs->x[j],rho) ;

		// Evaluate
		printf(" Supplied: ");
		for (int k = 0; k < gVarSize; k++)
			printf(" %.4f", xs->x[k]);
		printf("\n");

		xs->f = gQualityF(xs->x);
		fnew[count - 1] = xs->f;
		
		// get the worst score, and the average
		fworstit[count - 1] = 0.0;
		favg[count - 1] = 0.0;
		for (int i = 0; i < gPopSize; i++)
		{
			favg[count - 1] += gPopulation[i].f;
			if (gPopulation[i].f > fworstit[count - 1])
				fworstit[count - 1] = gPopulation[i].f;
		}
		favg[count - 1] /= gPopSize;
		
		// plot the current scores in the debug graph
		DeleteGraphPlot(debugpanel, debuggraph1, -1, VAL_DELAYED_DRAW);
		PlotY(debugpanel, debuggraph1, fnew,    count, VAL_DOUBLE, VAL_SCATTER,   VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0xFF0000);
		PlotY(debugpanel, debuggraph1, fbestit, count, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0x00FF00);
		PlotY(debugpanel, debuggraph1, fworstit, count, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0x0000FF);
		PlotY(debugpanel, debuggraph1, favg,    count, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0xFFFFFF);			
		//PlotY(debugpanel, debuggraph1, xs->x, gVarSize, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0x00FF00);			
		RefreshGraph(debugpanel, debuggraph1);
		
		/*printf("new genes:");
		for (int p=0; p < gVarSize; p++)
			printf(" %.3f", xs->x[p]);
		printf("\n");*/

		// See if it is better than the selected candidate for replacement
		if( xs->f < xi->f )
		{
			// yes, replace the selected 'bad' candidate
			for(unsigned int j=0 ; j<gVarSize ; j++) 
				xi->x[j] = xs->x[j];
			xi->f = xs->f ;

			// see if it is even better than the best specimen so far
			if( xi->f < fbest )
			{
				fbest = xi->f ;
				ibest = gPopulation-xi ;
				printf( "Best %d: %f\n", count, fbest ) ;
			}
		}
	}

	free(fnew);
	free(fbestit);

	return 0;
}


void GeneticDeinitialise()
{
	// Cleanup
	for(unsigned int i=0 ; i<gPopSize+1 ; i++) 
		free(gPopulation[i].x) ;
	free(gPopulation);
}
