//==============================================================================
//
// Title:       Genetic Algorithm header file
// Purpose:     Algorithms for performing genetic optimization
//
// Created on:  14-09-2012 by Rick van Bijnen, Bas van der Geer
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================


// variable type of the 'genes'
typedef double tGeneticVar ; 

// Pointer to function type
typedef double(*dFpd)(double*);

// initialise the population
void GeneticInitialize(int popsize, unsigned int varsize, tGeneticVar maxvariance, dFpd f);

// perform genetic optimization
int GeneticOptimize(int NumEval, int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6, int debuggraph1);

// cleanup
void GeneticDeinitialise();


