/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


void computeAderIntegration() {
  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(m_lts.dofs);
  real**                buffers                       = layer.var(m_lts.buffers);
  real**                derivatives                   = layer.var(m_lts.derivatives);
  LocalIntegrationData* localIntegration              = layer.var(m_lts.localIntegration);
  
#ifdef _OPENMP
  #pragma omp parallel 
  {
  #pragma omp for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    m_timeKernel.computeAder(              m_timeStepWidthSimulation,
                                           &m_globalData,
                                           &localIntegration[l_cell],
                                           dofs[l_cell],
                                           buffers[l_cell],
                                           derivatives[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeLocalWithoutAderIntegration() {
  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(m_lts.dofs);
  real**                buffers                       = layer.var(m_lts.buffers);
  LocalIntegrationData* localIntegration              = layer.var(m_lts.localIntegration);
  CellLocalInformation* cellInformation               = layer.var(m_lts.cellInformation);

#ifdef _OPENMP
  #pragma omp parallel 
  {
  #pragma omp for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    m_localKernel.computeIntegral(  cellInformation[l_cell].faceTypes,
                                    &m_globalData,
                                    &localIntegration[l_cell],
                                    buffers[l_cell],
                                    dofs[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeLocalIntegration() {
  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(m_lts.dofs);
  real**                buffers                       = layer.var(m_lts.buffers);
  real**                derivatives                   = layer.var(m_lts.derivatives);
  LocalIntegrationData* localIntegration              = layer.var(m_lts.localIntegration);
  CellLocalInformation* cellInformation               = layer.var(m_lts.cellInformation);

#ifdef _OPENMP
  #pragma omp parallel
  {
  #pragma omp for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    m_timeKernel.computeAder(      (double)m_timeStepWidthSimulation,
                                           &m_globalData,
                                           &localIntegration[l_cell],
                                           dofs[l_cell],
                                           buffers[l_cell],
                                           derivatives[l_cell] );

    m_localKernel.computeIntegral(        cellInformation[l_cell].faceTypes,
                                          &m_globalData,
                                          &localIntegration[l_cell],
                                           buffers[l_cell],
                                           dofs[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeNeighboringIntegration() {
  auto&                     layer                           = m_ltsTree.child(0).child<Interior>();
  unsigned                  nrOfCells                       = layer.getNumberOfCells();
  real                      (*dofs)[NUMBER_OF_ALIGNED_DOFS] = layer.var(m_lts.dofs);
  real*                     (*faceNeighbors)[4]             = layer.var(m_lts.faceNeighbors);
  CellDRMapping             (*drMapping)[4]                 = layer.var(m_lts.drMapping);
  NeighboringIntegrationData* neighboringIntegration        = layer.var(m_lts.neighboringIntegration);
  CellLocalInformation*       cellInformation               = layer.var(m_lts.cellInformation);
  
  real *l_timeIntegrated[4];
#ifdef ENABLE_MATRIX_PREFETCH
  real *l_faceNeighbors_prefetch[4];
#endif

#ifdef _OPENMP
#  ifdef ENABLE_MATRIX_PREFETCH
  #pragma omp parallel private(l_timeIntegrated, l_faceNeighbors_prefetch)
#  else
  #pragma omp parallel private(l_timeIntegrated)
#  endif
  {
  #pragma omp for schedule(static)
#endif
  for( int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    seissol::kernels::TimeCommon::computeIntegrals( m_timeKernel,
                                                    cellInformation[l_cell].ltsSetup,
                                                    cellInformation[l_cell].faceTypes,
                                                    0.0,
                                            (double)m_timeStepWidthSimulation,
                                                    faceNeighbors[l_cell],
#ifdef _OPENMP
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(&(m_globalData.integrationBufferLTS[omp_get_thread_num()*4*NUMBER_OF_ALIGNED_DOFS])),
#else
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(m_globalData.integrationBufferLTS),
#endif
                                                    l_timeIntegrated );

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != dynamicRupture) ? faceNeighbors[l_cell][1] : drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != dynamicRupture) ? faceNeighbors[l_cell][2] : drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != dynamicRupture) ? faceNeighbors[l_cell][3] : drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (nrOfCells-1) ) {
      l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != dynamicRupture) ? faceNeighbors[l_cell+1][0] : drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( cellInformation[l_cell].faceTypes,
                                               cellInformation[l_cell].faceRelations,
                                               drMapping[l_cell],
                                               &m_globalData,
                                               &neighboringIntegration[l_cell],
                                               l_timeIntegrated,
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_faceNeighbors_prefetch,
#endif
                                               dofs[l_cell]);
  }

#ifdef _OPENMP
  }
#endif
}

void computeDynRupGodunovState()
{
  seissol::initializers::Layer& layerData = m_dynRupTree.child(0).child<Interior>();
  DRFaceInformation*                    faceInformation                                                   = layerData.var(m_dynRup.faceInformation);
  DRGodunovData*                        godunovData                                                       = layerData.var(m_dynRup.godunovData);
  real**                                timeDerivativePlus                                                = layerData.var(m_dynRup.timeDerivativePlus);
  real**                                timeDerivativeMinus                                               = layerData.var(m_dynRup.timeDerivativeMinus);
  real                                (*godunov)[CONVERGENCE_ORDER][seissol::model::godunovState::reals]  = layerData.var(m_dynRup.godunov);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynRupKernel.computeGodunovState( faceInformation[face],
                                        &m_globalData,
                                       &godunovData[face],
                                        timeDerivativePlus[face],
                                        timeDerivativeMinus[face],
                                        godunov[face],
                                        timeDerivativePlus[prefetchFace],
                                        timeDerivativeMinus[prefetchFace] );
  }
}

