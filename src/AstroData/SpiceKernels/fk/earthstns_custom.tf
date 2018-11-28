KPL/FK
 
   FILE: earthstns_custom.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2018-11-28T11:02:12
   PINPOINT DEFINITIONS FILE: pinpoint_ex1.defs
   PINPOINT PCK FILE:         pck00010.tpc
   PINPOINT SPK FILE:         earthstns_custom.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'GROTTAGLIE'
   NAIF_BODY_CODE                      += 399109991
 
   NAIF_BODY_NAME                      += 'ROMA'
   NAIF_BODY_CODE                      += 399109992
 
   NAIF_BODY_NAME                      += 'MATERA'
   NAIF_BODY_CODE                      += 399109993
 
   NAIF_BODY_NAME                      += 'KIRUNA'
   NAIF_BODY_CODE                      += 399109994
 
   NAIF_BODY_NAME                      += 'TROLL'
   NAIF_BODY_CODE                      += 399109995
 
   NAIF_BODY_NAME                      += 'SVALDARD'
   NAIF_BODY_CODE                      += 399109996
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame GROTTAGLIE_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame GROTTAGLIE_TOPO is centered at the
      site GROTTAGLIE, which has Cartesian coordinates
 
         X (km):                  0.4631499414460E+04
         Y (km):                  0.1454383594070E+04
         Z (km):                  0.4123219994888E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        17.4333333333000
         Latitude  (deg):        40.5333333333000
         Altitude   (km):         0.7147000000213E-01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_GROTTAGLIE_TOPO               =  400109991
   FRAME_400109991_NAME                =  'GROTTAGLIE_TOPO'
   FRAME_400109991_CLASS               =  4
   FRAME_400109991_CLASS_ID            =  400109991
   FRAME_400109991_CENTER              =  399109991
 
   OBJECT_399109991_FRAME              =  'GROTTAGLIE_TOPO'
 
   TKFRAME_400109991_RELATIVE          =  'ITRF93'
   TKFRAME_400109991_SPEC              =  'ANGLES'
   TKFRAME_400109991_UNITS             =  'DEGREES'
   TKFRAME_400109991_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109991_ANGLES            =  (  -17.4333333333000,
                                             -49.4666666667000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame ROMA_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame ROMA_TOPO is centered at the
      site ROMA, which has Cartesian coordinates
 
         X (km):                  0.4642526773181E+04
         Y (km):                  0.1027759466320E+04
         Z (km):                  0.4236759467972E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        12.4827780000000
         Latitude  (deg):        41.8930560000000
         Altitude   (km):        -0.1380799999751E-01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_ROMA_TOPO                     =  400109992
   FRAME_400109992_NAME                =  'ROMA_TOPO'
   FRAME_400109992_CLASS               =  4
   FRAME_400109992_CLASS_ID            =  400109992
   FRAME_400109992_CENTER              =  399109992
 
   OBJECT_399109992_FRAME              =  'ROMA_TOPO'
 
   TKFRAME_400109992_RELATIVE          =  'ITRF93'
   TKFRAME_400109992_SPEC              =  'ANGLES'
   TKFRAME_400109992_UNITS             =  'DEGREES'
   TKFRAME_400109992_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109992_ANGLES            =  (  -12.4827780000000,
                                             -48.1069440000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame MATERA_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame MATERA_TOPO is centered at the
      site MATERA, which has Cartesian coordinates
 
         X (km):                  0.4641878479297E+04
         Y (km):                  0.1393027697827E+04
         Z (km):                  0.4133220167462E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        16.7045002000000
         Latitude  (deg):        40.6491012600000
         Altitude   (km):         0.4367000000010E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_MATERA_TOPO                   =  400109993
   FRAME_400109993_NAME                =  'MATERA_TOPO'
   FRAME_400109993_CLASS               =  4
   FRAME_400109993_CLASS_ID            =  400109993
   FRAME_400109993_CENTER              =  399109993
 
   OBJECT_399109993_FRAME              =  'MATERA_TOPO'
 
   TKFRAME_400109993_RELATIVE          =  'ITRF93'
   TKFRAME_400109993_SPEC              =  'ANGLES'
   TKFRAME_400109993_UNITS             =  'DEGREES'
   TKFRAME_400109993_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109993_ANGLES            =  (  -16.7045002000000,
                                             -49.3508987400000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame KIRUNA_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame KIRUNA_TOPO is centered at the
      site KIRUNA, which has Cartesian coordinates
 
         X (km):                  0.2251459728027E+04
         Y (km):                  0.8627615111859E+03
         Z (km):                  0.5885835596030E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        20.9668800000000
         Latitude  (deg):        67.8584280000000
         Altitude   (km):         0.7299999999995E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_KIRUNA_TOPO                   =  400109994
   FRAME_400109994_NAME                =  'KIRUNA_TOPO'
   FRAME_400109994_CLASS               =  4
   FRAME_400109994_CLASS_ID            =  400109994
   FRAME_400109994_CENTER              =  399109994
 
   OBJECT_399109994_FRAME              =  'KIRUNA_TOPO'
 
   TKFRAME_400109994_RELATIVE          =  'ITRF93'
   TKFRAME_400109994_SPEC              =  'ANGLES'
   TKFRAME_400109994_UNITS             =  'DEGREES'
   TKFRAME_400109994_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109994_ANGLES            =  (  -20.9668800000000,
                                             -22.1415720000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame TROLL_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame TROLL_TOPO is centered at the
      site TROLL, which has Cartesian coordinates
 
         X (km):                  0.1973469690474E+04
         Y (km):                  0.8731388348027E+02
         Z (km):                 -0.6044956617609E+04
 
      and planetodetic coordinates
 
         Longitude (deg):         2.5333333300000
         Latitude  (deg):       -72.0167000000000
         Altitude   (km):         0.7308529999996E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_TROLL_TOPO                    =  400109995
   FRAME_400109995_NAME                =  'TROLL_TOPO'
   FRAME_400109995_CLASS               =  4
   FRAME_400109995_CLASS_ID            =  400109995
   FRAME_400109995_CENTER              =  399109995
 
   OBJECT_399109995_FRAME              =  'TROLL_TOPO'
 
   TKFRAME_400109995_RELATIVE          =  'ITRF93'
   TKFRAME_400109995_SPEC              =  'ANGLES'
   TKFRAME_400109995_UNITS             =  'DEGREES'
   TKFRAME_400109995_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109995_ANGLES            =  (   -2.5333333300000,
                                            -162.0167000000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame SVALDARD_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame SVALDARD_TOPO is centered at the
      site SVALDARD, which has Cartesian coordinates
 
         X (km):                  0.1258485308372E+04
         Y (km):                  0.3468286061818E+03
         Z (km):                  0.6222934840953E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        15.4077860000000
         Latitude  (deg):        78.2297720000000
         Altitude   (km):         0.7308530000022E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_SVALDARD_TOPO                 =  400109996
   FRAME_400109996_NAME                =  'SVALDARD_TOPO'
   FRAME_400109996_CLASS               =  4
   FRAME_400109996_CLASS_ID            =  400109996
   FRAME_400109996_CENTER              =  399109996
 
   OBJECT_399109996_FRAME              =  'SVALDARD_TOPO'
 
   TKFRAME_400109996_RELATIVE          =  'ITRF93'
   TKFRAME_400109996_SPEC              =  'ANGLES'
   TKFRAME_400109996_UNITS             =  'DEGREES'
   TKFRAME_400109996_AXES              =  ( 3, 2, 3 )
   TKFRAME_400109996_ANGLES            =  (  -15.4077860000000,
                                             -11.7702280000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file pinpoint_ex1.defs
--------------------------------------------------------------------------------
 
 
 
Start of PINPOINT inputs.
 
 
begindata
 
 
   SITES                   += ('GROTTAGLIE','ROMA', 'MATERA','KIRUNA','TROLL','SVALDARD')
 
 
   GROTTAGLIE_CENTER        = 399
   GROTTAGLIE_FRAME         = 'ITRF93'
   GROTTAGLIE_IDCODE        = 399109991
   GROTTAGLIE_LATLON        = (  40.5333333333,  17.4333333333,  0.071470 )
   GROTTAGLIE_UP            = 'Z'
   GROTTAGLIE_NORTH         = 'X'
 
 
   ROMA_CENTER              = 399
   ROMA_FRAME               = 'ITRF93'
   ROMA_IDCODE              = 399109992
   ROMA_LATLON              = (   41.893056, 12.482778, -0.013808 )
   ROMA_UP                  = 'Z'
   ROMA_NORTH               = 'X'
 
 
   MATERA_CENTER            = 399
   MATERA_FRAME             = 'ITRF93'
   MATERA_IDCODE            = 399109993
   MATERA_LATLON            = (  40.64910126, 16.70450020,  0.436700 )
   MATERA_UP                = 'Z'
   MATERA_NORTH             = 'X'
 
 
   KIRUNA_CENTER            = 399
   KIRUNA_FRAME             = 'ITRF93'
   KIRUNA_IDCODE            = 399109994
   KIRUNA_LATLON            = ( 67.858428, 20.966880,  0.730000 )
   KIRUNA_UP                = 'Z'
   KIRUNA_NORTH             = 'X'
 
 
   TROLL_CENTER             = 399
   TROLL_FRAME              = 'ITRF93'
   TROLL_IDCODE             = 399109995
   TROLL_LATLON             = ( -72.0167, 2.53333333,  0.730853 )
   TROLL_UP                 = 'Z'
   TROLL_NORTH              = 'X'
 
 
   SVALDARD_CENTER          = 399
   SVALDARD_FRAME           = 'ITRF93'
   SVALDARD_IDCODE          = 399109996
   SVALDARD_LATLON          = ( 78.229772, 15.407786,  0.730853 )
   SVALDARD_UP              = 'Z'
   SVALDARD_NORTH           = 'X'
 
 
begintext
 
 
 
End of PINPOINT inputs.
 
 
begintext
 
[End of definitions file]
 
