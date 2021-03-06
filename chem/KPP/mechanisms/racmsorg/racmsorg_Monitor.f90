! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Utility Data Module File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : racmsorg_Monitor.f90
! Time                 : Tue Jan 11 14:33:12 2022
! Working directory    : /jathar-scratch/yicongh/WRF_GoAmazon_gasaerochem.DEV/chem/KPP/mechanisms/racmsorg
! Equation file        : racmsorg.kpp
! Output root filename : racmsorg
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racmsorg_Monitor


  CHARACTER(LEN=12), PARAMETER, DIMENSION(75) :: SPC_NAMES = (/ &
     'SULF        ','CO2         ','ORA1        ', &
     'ORA2        ','SO2         ','O1D         ', &
     'HC5         ','TOL         ','XYL         ', &
     'N2O5        ','HC8         ','HC3         ', &
     'ETH         ','CH4         ','UDD         ', &
     'HNO4        ','OP1         ','HONO        ', &
     'H2O2        ','PHO         ','ADDT        ', &
     'ADDX        ','HKET        ','ETE         ', &
     'ADDC        ','PAA         ','HNO3        ', &
     'CO          ','API         ','LIM         ', &
     'PAN         ','CSL         ','DIEN        ', &
     'GLY         ','TPAN        ','ETEP        ', &
     'ISO         ','OLTP        ','OLIP        ', &
     'MGLY        ','CSLP        ','KET         ', &
     'LIMP        ','HC5P        ','HC8P        ', &
     'TOLP        ','XYLP        ','HCHO        ', &
     'APIP        ','ISOP        ','MACR        ', &
     'HC3P        ','ALD         ','DCB         ', &
     'TCO3        ','XO2         ','OLT         ', &
     'OLI         ','OLNN        ','OLND        ', &
     'ETHP        ','O3P         ','O3          ', &
     'KETP        ','MO2         ','ACO3        ', &
     'HO2         ','HO          ','ONIT        ', &
     'NO3         ','OP2         ','NO          ', &
     'NO2         ','H2O         ','M           ' /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=12), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_0 = (/ &
     '        NO2 --> O3P + NO                                                                            ', &
     '         O3 --> O1D                                                                                 ', &
     '         O3 --> O3P                                                                                 ', &
     '       HONO --> HO + NO                                                                             ', &
     '       HNO3 --> HO + NO2                                                                            ', &
     '       HNO4 --> 0.65 HO2 + 0.35 HO + 0.35 NO3 + 0.65 NO2                                            ', &
     '        NO3 --> NO                                                                                  ', &
     '        NO3 --> O3P + NO2                                                                           ', &
     '       H2O2 --> 2 HO                                                                                ', &
     '       HCHO --> CO                                                                                  ', &
     '       HCHO --> CO + 2 HO2                                                                          ', &
     '        ALD --> CO + MO2 + HO2                                                                      ', &
     '        OP1 --> HCHO + HO2 + HO                                                                     ', &
     '        OP2 --> ALD + HO2 + HO                                                                      ', &
     '        PAA --> MO2 + HO                                                                            ', &
     '        KET --> ETHP + ACO3                                                                         ', &
     '        GLY --> 1.87 CO + 0.13 HCHO                                                                 ', &
     '        GLY --> 1.55 CO + 0.45 HCHO + 0.8 HO2                                                       ', &
     '       MGLY --> CO + ACO3 + HO2                                                                     ', &
     '        DCB --> TCO3 + HO2                                                                          ', &
     '       ONIT --> 0.8 KET + 0.2 ALD + HO2 + NO2                                                       ', &
     '       MACR --> CO + HCHO + ACO3 + HO2                                                              ', &
     '       HKET --> HCHO + ACO3 + HO2                                                                   ', &
     '    O3P + M --> O3                                                                                  ', &
     '   O3P + O3 --> M                                                                                   ', &
     '    O1D + M --> O3                                                                                  ', &
     '    O1D + M --> O3P                                                                                 ', &
     '  O1D + H2O --> 2 HO                                                                                ', &
     '    O3 + HO --> HO2                                                                                 ', &
     '   O3 + HO2 --> HO                                                                                  ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_1 = (/ &
     '   HO2 + HO --> H2O                                                                                 ', &
     '  H2O2 + HO --> HO2 + H2O                                                                           ', &
     '      2 HO2 --> H2O2                                                                                ', &
     '2 HO2 + H2O --> H2O2 + H2O                                                                          ', &
     '   O3P + NO --> NO2                                                                                 ', &
     '  O3P + NO2 --> NO                                                                                  ', &
     '  O3P + NO2 --> NO3                                                                                 ', &
     '    HO + NO --> HONO                                                                                ', &
     '   HO + NO2 --> HNO3                                                                                ', &
     '   HO + NO3 --> HO2 + NO2                                                                           ', &
     '   HO2 + NO --> HO + NO2                                                                            ', &
     '  HO2 + NO2 --> HNO4                                                                                ', &
     '       HNO4 --> HO2 + NO2                                                                           ', &
     '  HO2 + NO3 --> 0.3 HNO3 + 0.7 HO + 0.7 NO2                                                         ', &
     '  HONO + HO --> NO2 + H2O                                                                           ', &
     '  HNO3 + HO --> NO3 + H2O                                                                           ', &
     '  HNO4 + HO --> NO2 + H2O                                                                           ', &
     '    O3 + NO --> NO2                                                                                 ', &
     '   O3 + NO2 --> NO3                                                                                 ', &
     '   2 NO + M --> 2 NO2                                                                               ', &
     '   NO3 + NO --> 2 NO2                                                                               ', &
     '  NO3 + NO2 --> NO + NO2                                                                            ', &
     '  NO3 + NO2 --> N2O5                                                                                ', &
     '       N2O5 --> NO3 + NO2                                                                           ', &
     '      2 NO3 --> 2 NO2                                                                               ', &
     '     HO + M --> HO2 + H2O                                                                           ', &
     '   SO2 + HO --> SULF + HO2                                                                          ', &
     '    CO + HO --> CO2 + HO2                                                                           ', &
     '  ISO + O3P --> 0.01 CO + 0.05 HCHO + 0.13 DCB + 0.15 XO2 + 0.86 OLT + 0.28 HO2 + 0.02 HO           ', &
     ' MACR + O3P --> ALD                                                                                 ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_2 = (/ &
     '   CH4 + HO --> MO2 + H2O                                                                           ', &
     '   ETH + HO --> ETHP + H2O                                                                          ', &
     '   HC3 + HO --> 0.036 ORA1 + 0.036 CO + 0.036 GLY + 0.01 HCHO + 0.583 HC3P + 0.335 ALD + 0.381 HO2 +', &
     '   HC5 + HO --> 0.25 KET + 0.75 HC5P + 0.25 HO2 + H2O                                               ', &
     '   HC8 + HO --> 0.024 HKET + 0.9511 HC8P + 0.025 ALD + 0.049 HO2 + H2O                              ', &
     '   ETE + HO --> ETEP                                                                                ', &
     '   OLT + HO --> OLTP                                                                                ', &
     '   OLI + HO --> OLIP                                                                                ', &
     '  DIEN + HO --> ISOP                                                                                ', &
     '   ISO + HO --> ISOP                                                                                ', &
     '   API + HO --> APIP                                                                                ', &
     '   LIM + HO --> LIMP                                                                                ', &
     '   TOL + HO --> 0.9 ADDT + 0.1 XO2 + 0.1 HO2                                                        ', &
     '   XYL + HO --> 0.9 ADDX + 0.1 XO2 + 0.1 HO2                                                        ', &
     '   CSL + HO --> 0.1 PHO + 0.85 ADDC + 0.05 XO2 + 0.05 HO2                                           ', &
     '  HCHO + HO --> CO + HO2 + H2O                                                                      ', &
     '   ALD + HO --> ACO3 + H2O                                                                          ', &
     '   KET + HO --> KETP + H2O                                                                          ', &
     '  HKET + HO --> MGLY + HO2 + H2O                                                                    ', &
     '   GLY + HO --> 2 CO + HO2 + H2O                                                                    ', &
     '  MGLY + HO --> CO + ACO3 + H2O                                                                     ', &
     '  MACR + HO --> 0.41 HKET + 0.41 CO + 0.08 MGLY + 0.08 HCHO + 0.51 TCO3 + 0.49 XO2 + 0.49 HO2       ', &
     '   DCB + HO --> 0.35 UDD + 0.15 GLY + 0.15 MGLY + 0.5 TCO3 + 0.5 XO2 + 0.5 HO2                      ', &
     '   UDD + HO --> 0.12 KET + 0.88 ALD + HO2                                                           ', &
     '   OP1 + HO --> 0.35 HCHO + 0.65 MO2 + 0.35 HO                                                      ', &
     '   HO + OP2 --> 0.41 KET + 0.44 HC3P + 0.08 ALD + 0.07 XO2 + 0.49 HO                                ', &
     '   PAA + HO --> 0.35 HCHO + 0.35 XO2 + 0.65 ACO3 + 0.35 HO2                                         ', &
     '   PAN + HO --> HCHO + XO2 + NO3 + H2O                                                              ', &
     '  TPAN + HO --> 0.6 HKET + 0.4 PAN + 0.4 HCHO + XO2 + 0.4 HO2 + 0.6 NO3                             ', &
     '  HO + ONIT --> HC3P + NO2 + H2O                                                                    ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_3 = (/ &
     ' HCHO + NO3 --> HNO3 + CO + HO2                                                                     ', &
     '  ALD + NO3 --> HNO3 + ACO3                                                                         ', &
     '  GLY + NO3 --> HNO3 + 2 CO + HO2                                                                   ', &
     ' MGLY + NO3 --> HNO3 + CO + ACO3                                                                    ', &
     ' MACR + NO3 --> 0.2 HNO3 + 0.8 CO + 0.2 TCO3 + 0.8 OLNN                                             ', &
     '  DCB + NO3 --> 0.5 HNO3 + 0.25 GLY + 0.25 MGLY + 0.03 KET + 0.25 ALD + 0.5 TCO3 + 0.5 XO2 + 0.5 HO2', &
     '  CSL + NO3 --> PHO + HNO3                                                                          ', &
     '  ETE + NO3 --> 0.8 OLNN + 0.2 OLND                                                                 ', &
     '  OLT + NO3 --> 0.43 OLNN + 0.57 OLND                                                               ', &
     '  OLI + NO3 --> 0.11 OLNN + 0.89 OLND                                                               ', &
     ' DIEN + NO3 --> 0.9 MACR + 0.9 OLNN + 0.1 OLND                                                      ', &
     '  ISO + NO3 --> 0.9 MACR + 0.9 OLNN + 0.1 OLND                                                      ', &
     '  API + NO3 --> 0.1 OLNN + 0.9 OLND                                                                 ', &
     '  LIM + NO3 --> 0.13 OLNN + 0.87 OLND                                                               ', &
     ' TPAN + NO3 --> 0.4 PAN + 0.4 HCHO + XO2 + 0.6 ONIT + 0.6 NO3 + 0.4 NO2                             ', &
     '   ETE + O3 --> 0.37 ORA1 + 0.43 CO + HCHO + 0.26 HO2 + 0.12 HO                                     ', &
     '   OLT + O3 --> 0.14 ORA1 + 0.1 ORA2 + 0.03 ETH + 0.06 CH4 + 0.006 H2O2 + 0.37 CO + 0.03 KET + 0.64 ', &
     '   OLI + O3 --> 0.14 ORA2 + 0.06 ETH + 0.07 CH4 + 0.011 H2O2 + 0.3 CO + 0.16 KET + 0.02 HCHO + 0.99 ', &
     '  DIEN + O3 --> 0.15 ORA1 + 0.001 H2O2 + 0.36 CO + 0.9 HCHO + 0.39 MACR + 0.13 XO2 + 0.35 OLT + 0.09', &
     '   ISO + O3 --> 0.15 ORA1 + 0.001 H2O2 + 0.36 CO + 0.9 HCHO + 0.39 MACR + 0.13 XO2 + 0.35 OLT + 0.09', &
     '   API + O3 --> 0.02 H2O2 + 0.14 CO + 0.53 KET + 0.65 ALD + 0.2 ETHP + 0.42 KETP + 0.1 HO2 + 0.85 HO', &
     '   LIM + O3 --> 0.01 ORA1 + 0.07 ORA2 + 0.02 H2O2 + 0.14 CO + 0.04 HCHO + 0.79 MACR + 0.46 OLT + 0.1', &
     '  MACR + O3 --> 0.22 ORA1 + 0.13 ORA2 + 0.54 CO + 0.6 MGLY + 0.4 HCHO + 0.13 ACO3 + 0.29 HO2 + 0.07 ', &
     '   DCB + O3 --> 0.11 ORA1 + 0.21 ORA2 + 0.11 PAA + 0.66 CO + 0.5 GLY + 0.62 MGLY + 0.16 ALD + 0.28 A', &
     '  TPAN + O3 --> 0.11 ORA1 + 0.13 CO + 0.3 PAN + 0.7 HCHO + 0.7 ACO3 + 0.08 HO2 + 0.036 HO + 0.7 NO2 ', &
     '  PHO + NO2 --> 0.1 CSL + ONIT                                                                      ', &
     '  PHO + HO2 --> CSL                                                                                 ', &
     ' ADDT + NO2 --> HONO + CSL                                                                          ', &
     '   ADDT + M --> 0.02 CSL + 0.98 TOLP + 0.02 HO2                                                     ', &
     '  ADDT + O3 --> CSL + HO                                                                            ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_4 = (/ &
     ' ADDX + NO2 --> HONO + CSL                                                                          ', &
     '   ADDX + M --> 0.02 CSL + 0.98 XYLP + 0.02 HO2                                                     ', &
     '  ADDX + O3 --> CSL + HO                                                                            ', &
     ' ADDC + NO2 --> HONO + CSL                                                                          ', &
     '   ADDC + M --> 0.02 CSL + 0.98 CSLP + 0.02 HO2                                                     ', &
     '  ADDC + O3 --> CSL + HO                                                                            ', &
     ' ACO3 + NO2 --> PAN                                                                                 ', &
     '        PAN --> ACO3 + NO2                                                                          ', &
     ' TCO3 + NO2 --> TPAN                                                                                ', &
     '       TPAN --> TCO3 + NO2                                                                          ', &
     '   MO2 + NO --> HCHO + HO2 + NO2                                                                    ', &
     '  ETHP + NO --> ALD + HO2 + NO2                                                                     ', &
     '  HC3P + NO --> 0.063 GLY + 0.623 KET + 0.047 HCHO + 0.233 ALD + 0.048 XO2 + 0.048 ETHP + 0.15 MO2 +', &
     '  HC5P + NO --> 0.722 KET + 0.021 HCHO + 0.211 ALD + 0.334 XO2 + 0.245 ETHP + 0.031 MO2 + 0.599 HO2 ', &
     '  HC8P + NO --> 0.642 KET + 0.15 ALD + 0.416 XO2 + 0.133 ETHP + 0.606 HO2 + 0.261 ONIT + 0.739 NO2  ', &
     '  ETEP + NO --> 1.6 HCHO + 0.2 ALD + HO2 + NO2                                                      ', &
     '  OLTP + NO --> 0.06 KET + HCHO + 0.94 ALD + HO2 + NO2                                              ', &
     '  OLIP + NO --> 0.29 KET + 1.71 ALD + HO2 + NO2                                                     ', &
     '  ISOP + NO --> 0.606 HCHO + 0.446 MACR + 0.354 OLT + 0.847 HO2 + 0.153 ONIT + 0.847 NO2            ', &
     '  APIP + NO --> 0.8 KET + 0.8 ALD + 0.8 HO2 + 0.2 ONIT + 0.8 NO2                                    ', &
     '  LIMP + NO --> 0.25 HCHO + 0.4 MACR + 0.25 OLI + 0.65 HO2 + 0.35 ONIT + 0.65 NO2                   ', &
     '  TOLP + NO --> 1.2 GLY + 0.65 MGLY + 0.5 DCB + 0.95 HO2 + 0.05 ONIT + 0.95 NO2                     ', &
     '  XYLP + NO --> 0.35 GLY + 0.6 MGLY + 0.95 DCB + 0.95 HO2 + 0.05 ONIT + 0.95 NO2                    ', &
     '  CSLP + NO --> GLY + MGLY + HO2 + NO2                                                              ', &
     '  ACO3 + NO --> MO2 + NO2                                                                           ', &
     '  TCO3 + NO --> HCHO + ACO3 + NO2                                                                   ', &
     '  KETP + NO --> 0.54 MGLY + 0.46 ALD + 0.16 XO2 + 0.23 ACO3 + 0.77 HO2 + NO2                        ', &
     '  OLNN + NO --> HO2 + ONIT + NO2                                                                    ', &
     '  OLND + NO --> 0.464 KET + 0.287 HCHO + 1.24 ALD + 2 NO2                                           ', &
     '  MO2 + HO2 --> OP1                                                                                 ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_5 = (/ &
     ' ETHP + HO2 --> OP2                                                                                 ', &
     ' HC3P + HO2 --> OP2                                                                                 ', &
     ' HC5P + HO2 --> OP2                                                                                 ', &
     ' HC8P + HO2 --> OP2                                                                                 ', &
     ' ETEP + HO2 --> OP2                                                                                 ', &
     ' OLTP + HO2 --> OP2                                                                                 ', &
     ' OLIP + HO2 --> OP2                                                                                 ', &
     ' ISOP + HO2 --> OP2                                                                                 ', &
     ' APIP + HO2 --> OP2                                                                                 ', &
     ' LIMP + HO2 --> OP2                                                                                 ', &
     ' TOLP + HO2 --> OP2                                                                                 ', &
     ' XYLP + HO2 --> OP2                                                                                 ', &
     ' CSLP + HO2 --> OP2                                                                                 ', &
     ' ACO3 + HO2 --> PAA                                                                                 ', &
     ' ACO3 + HO2 --> ORA2 + O3                                                                           ', &
     ' TCO3 + HO2 --> OP2                                                                                 ', &
     ' TCO3 + HO2 --> ORA2 + O3                                                                           ', &
     ' KETP + HO2 --> OP2                                                                                 ', &
     ' OLNN + HO2 --> ONIT                                                                                ', &
     ' OLND + HO2 --> ONIT                                                                                ', &
     '      2 MO2 --> 1.33 HCHO + 0.66 HO2                                                                ', &
     ' ETHP + MO2 --> 0.75 HCHO + 0.75 ALD + HO2                                                          ', &
     ' HC3P + MO2 --> 0.119 GLY + 0.005 MGLY + 0.018 KET + 0.81 HCHO + 0.58 ALD + 0.085 XO2 + 0.007 MO2 + ', &
     ' HC5P + MO2 --> 0.24 KET + 0.829 HCHO + 0.523 ALD + 0.245 XO2 + 0.014 ETHP + 0.049 MO2 + 0.946 HO2  ', &
     ' HC8P + MO2 --> 0.419 KET + 0.753 HCHO + 0.411 ALD + 0.322 XO2 + 0.013 ETHP + 0.993 HO2             ', &
     ' ETEP + MO2 --> 1.55 HCHO + 0.35 ALD + HO2                                                          ', &
     ' OLTP + MO2 --> 0.081 KET + 1.25 HCHO + 0.669 ALD + HO2                                             ', &
     ' OLIP + MO2 --> 0.313 KET + 0.755 HCHO + 0.932 ALD + HO2                                            ', &
     ' ISOP + MO2 --> 1.09 HCHO + 0.55 MACR + 0.37 OLT + 0.08 OLI + HO2                                   ', &
     ' APIP + MO2 --> KET + HCHO + ALD + 2 HO2                                                            ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_6 = (/ &
     ' LIMP + MO2 --> 1.4 HCHO + 0.6 MACR + 0.4 OLI + 2 HO2                                               ', &
     ' TOLP + MO2 --> 0.65 GLY + 0.35 MGLY + HCHO + DCB + HO2                                             ', &
     ' XYLP + MO2 --> 0.37 GLY + 0.63 MGLY + HCHO + DCB + HO2                                             ', &
     ' CSLP + MO2 --> GLY + MGLY + HCHO + 2 HO2                                                           ', &
     ' MO2 + ACO3 --> HCHO + MO2 + HO2                                                                    ', &
     ' MO2 + ACO3 --> ORA2 + HCHO                                                                         ', &
     ' TCO3 + MO2 --> 2 HCHO + ACO3 + HO2                                                                 ', &
     ' TCO3 + MO2 --> ORA2 + HCHO                                                                         ', &
     ' KETP + MO2 --> 0.3 HKET + 0.4 MGLY + 0.75 HCHO + 0.3 ALD + 0.08 XO2 + 0.12 ACO3 + 0.88 HO2         ', &
     ' OLNN + MO2 --> 0.75 HCHO + HO2 + ONIT                                                              ', &
     ' OLND + MO2 --> 0.149 KET + 0.96 HCHO + 0.64 ALD + 0.5 HO2 + 0.5 ONIT + 0.5 NO2                     ', &
     'ETHP + ACO3 --> 0.5 ORA2 + ALD + 0.5 MO2 + 0.5 HO2                                                  ', &
     'HC3P + ACO3 --> 0.499 ORA2 + 0.1 GLY + 0.004 MGLY + 0.127 KET + 0.091 HCHO + 0.724 ALD + 0.071 XO2 +', &
     'HC5P + ACO3 --> 0.495 ORA2 + 0.33 KET + 0.076 HCHO + 0.677 ALD + 0.237 XO2 + 0.018 ETHP + 0.554 MO2 ', &
     'HC8P + ACO3 --> 0.495 ORA2 + 0.581 KET + 0.497 ALD + 0.318 XO2 + 0.015 ETHP + 0.507 MO2 + 0.489 HO2 ', &
     'ETEP + ACO3 --> 0.5 ORA2 + 0.8 HCHO + 0.6 ALD + 0.5 MO2 + 0.5 HO2                                   ', &
     'OLTP + ACO3 --> 0.499 ORA2 + 0.141 KET + 0.501 HCHO + 0.859 ALD + 0.501 MO2 + 0.501 HO2             ', &
     'OLIP + ACO3 --> 0.49 ORA2 + 0.569 KET + 0.941 ALD + 0.51 MO2 + 0.51 HO2                             ', &
     'ISOP + ACO3 --> 0.494 ORA2 + 0.34 HCHO + 0.771 MACR + 0.229 OLT + 0.506 MO2 + 0.506 HO2             ', &
     'APIP + ACO3 --> KET + ALD + MO2 + HO2                                                               ', &
     'LIMP + ACO3 --> 0.4 HCHO + 0.6 MACR + 0.4 OLI + MO2 + HO2                                           ', &
     'TOLP + ACO3 --> 0.65 GLY + 0.35 MGLY + DCB + MO2 + HO2                                              ', &
     'XYLP + ACO3 --> 0.37 GLY + 0.63 MGLY + DCB + MO2 + HO2                                              ', &
     'CSLP + ACO3 --> GLY + MGLY + MO2 + HO2                                                              ', &
     '     2 ACO3 --> 2 MO2                                                                               ', &
     'TCO3 + ACO3 --> HCHO + MO2 + ACO3                                                                   ', &
     'KETP + ACO3 --> 0.5 ORA2 + 0.54 MGLY + 0.11 KET + 0.35 ALD + 0.08 XO2 + 0.5 MO2 + 0.12 ACO3 + 0.38 H', &
     'OLNN + ACO3 --> 0.5 ORA2 + 0.5 MO2 + 0.5 HO2 + ONIT                                                 ', &
     'OLND + ACO3 --> 0.484 ORA2 + 0.167 KET + 0.207 HCHO + 0.65 ALD + 0.516 MO2 + 0.484 ONIT + 0.516 NO2 ', &
     '     2 OLNN --> HO2 + 2 ONIT                                                                        ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(27) :: EQN_NAMES_7 = (/ &
     'OLNN + OLND --> 0.149 KET + 0.202 HCHO + 0.64 ALD + 0.5 HO2 + 1.5 ONIT + 0.5 NO2                    ', &
     '     2 OLND --> 0.285 KET + 0.504 HCHO + 1.21 ALD + ONIT + NO2                                      ', &
     '  MO2 + NO3 --> HCHO + HO2 + NO2                                                                    ', &
     ' ETHP + NO3 --> ALD + HO2 + NO2                                                                     ', &
     ' HC3P + NO3 --> 0.063 GLY + 0.67 KET + 0.048 HCHO + 0.243 ALD + 0.051 XO2 + 0.053 ETHP + 0.155 MO2 +', &
     ' HC5P + NO3 --> 0.828 KET + 0.021 HCHO + 0.239 ALD + 0.391 XO2 + 0.262 ETHP + 0.04 MO2 + 0.699 HO2 +', &
     ' HC8P + NO3 --> 0.88 KET + 0.187 ALD + 0.587 XO2 + 0.155 ETHP + 0.845 HO2 + NO2                     ', &
     ' ETEP + NO3 --> 1.6 HCHO + 0.2 ALD + HO2 + NO2                                                      ', &
     ' OLTP + NO3 --> 0.06 KET + HCHO + 0.94 ALD + HO2 + NO2                                              ', &
     ' OLIP + NO3 --> 0.29 KET + 1.71 ALD + HO2 + NO2                                                     ', &
     ' ISOP + NO3 --> 0.686 HCHO + 0.6 MACR + 0.4 OLT + HO2 + NO2                                         ', &
     ' APIP + NO3 --> KET + ALD + HO2 + NO2                                                               ', &
     ' LIMP + NO3 --> 0.4 HCHO + 0.6 MACR + 0.4 OLI + HO2 + NO2                                           ', &
     ' TOLP + NO3 --> 1.3 GLY + 0.7 MGLY + 0.5 DCB + HO2 + NO2                                            ', &
     ' XYLP + NO3 --> 0.74 GLY + 1.26 MGLY + DCB + HO2 + NO2                                              ', &
     ' CSLP + NO3 --> GLY + MGLY + HO2 + NO2                                                              ', &
     ' ACO3 + NO3 --> MO2 + NO2                                                                           ', &
     ' TCO3 + NO3 --> HCHO + ACO3 + NO2                                                                   ', &
     ' KETP + NO3 --> 0.54 MGLY + 0.46 ALD + 0.16 XO2 + 0.23 ACO3 + 0.77 HO2 + NO2                        ', &
     ' OLNN + NO3 --> HO2 + ONIT + NO2                                                                    ', &
     ' OLND + NO3 --> 0.469 KET + 0.28 HCHO + 1.24 ALD + 2 NO2                                            ', &
     '  XO2 + HO2 --> OP2                                                                                 ', &
     '  XO2 + MO2 --> HCHO + HO2                                                                          ', &
     ' XO2 + ACO3 --> MO2                                                                                 ', &
     '      2 XO2 --> M                                                                                   ', &
     '   XO2 + NO --> NO2                                                                                 ', &
     '  XO2 + NO3 --> NO2                                                                                 ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(237) :: EQN_NAMES = (/&
    EQN_NAMES_0, EQN_NAMES_1, EQN_NAMES_2, EQN_NAMES_3, EQN_NAMES_4, &
    EQN_NAMES_5, EQN_NAMES_6, EQN_NAMES_7 /)

! INLINED global variables

! End INLINED global variables


END MODULE racmsorg_Monitor
