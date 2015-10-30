# Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
# retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.  
# 
#     * Neither the name of Sandia Corporation nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

# Microsoft Developer Studio Project File - Name="exodusii" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=exodusii - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "exodusii.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "exodusii.mak" CFG="exodusii - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "exodusii - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "exodusii - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "exodusii - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "$(CUBITROOT)/netcdf/netcdf-3.4.snl/include" /I "../../include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "EX_ERR_STR" /D "NT" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../lib/nt/libexoIIv2c406.lib"

!ELSEIF  "$(CFG)" == "exodusii - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "$(CUBITROOT)/netcdf/netcdf-3.4.snl/include" /I "../../include" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "EX_ERR_STR" /D "NT" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../lib/nt/libexoIIv2c406_db.lib"

!ENDIF 

# Begin Target

# Name "exodusii - Win32 Release"
# Name "exodusii - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\ex_conv.c
# End Source File
# Begin Source File

SOURCE=..\ex_utils.c
# End Source File
# Begin Source File

SOURCE=..\exclos.c
# End Source File
# Begin Source File

SOURCE=..\excn2s.c
# End Source File
# Begin Source File

SOURCE=..\excopy.c
# End Source File
# Begin Source File

SOURCE=..\excre.c
# End Source File
# Begin Source File

SOURCE=..\exerr.c
# End Source File
# Begin Source File

SOURCE=..\exgatm.c
# End Source File
# Begin Source File

SOURCE=..\exgcns.c
# End Source File
# Begin Source File

SOURCE=..\exgcon.c
# End Source File
# Begin Source File

SOURCE=..\exgcor.c
# End Source File
# Begin Source File

SOURCE=..\exgcss.c
# End Source File
# Begin Source File

SOURCE=..\exgeat.c
# End Source File
# Begin Source File

SOURCE=..\exgebi.c
# End Source File
# Begin Source File

SOURCE=..\exgelb.c
# End Source File
# Begin Source File

SOURCE=..\exgelc.c
# End Source File
# Begin Source File

SOURCE=..\exgem.c
# End Source File
# Begin Source File

SOURCE=..\exgenm.c
# End Source File
# Begin Source File

SOURCE=..\exgev.c
# End Source File
# Begin Source File

SOURCE=..\exgevt.c
# End Source File
# Begin Source File

SOURCE=..\exgfrm.c
# End Source File
# Begin Source File

SOURCE=..\exggv.c
# End Source File
# Begin Source File

SOURCE=..\exggvt.c
# End Source File
# Begin Source File

SOURCE=..\exginf.c
# End Source File
# Begin Source File

SOURCE=..\exgini.c
# End Source File
# Begin Source File

SOURCE=..\exgmap.c
# End Source File
# Begin Source File

SOURCE=..\exgmp.c
# End Source File
# Begin Source File

SOURCE=..\exgnm.c
# End Source File
# Begin Source File

SOURCE=..\exgnnm.c
# End Source File
# Begin Source File

SOURCE=..\exgnp.c
# End Source File
# Begin Source File

SOURCE=..\exgns.c
# End Source File
# Begin Source File

SOURCE=..\exgnsd.c
# End Source File
# Begin Source File

SOURCE=..\exgnsi.c
# End Source File
# Begin Source File

SOURCE=..\exgnv.c
# End Source File
# Begin Source File

SOURCE=..\exgnvt.c
# End Source File
# Begin Source File

SOURCE=..\exgp.c
# End Source File
# Begin Source File

SOURCE=..\exgpa.c
# End Source File
# Begin Source File

SOURCE=..\exgpn.c
# End Source File
# Begin Source File

SOURCE=..\exgqa.c
# End Source File
# Begin Source File

SOURCE=..\exgsnl.c
# End Source File
# Begin Source File

SOURCE=..\exgsp.c
# End Source File
# Begin Source File

SOURCE=..\exgss.c
# End Source File
# Begin Source File

SOURCE=..\exgssc.c
# End Source File
# Begin Source File

SOURCE=..\exgssd.c
# End Source File
# Begin Source File

SOURCE=..\exgssi.c
# End Source File
# Begin Source File

SOURCE=..\exgssn.c
# End Source File
# Begin Source File

SOURCE=..\exgtim.c
# End Source File
# Begin Source File

SOURCE=..\exgvan.c
# End Source File
# Begin Source File

SOURCE=..\exgvnm.c
# End Source File
# Begin Source File

SOURCE=..\exgvp.c
# End Source File
# Begin Source File

SOURCE=..\exgvtt.c
# End Source File
# Begin Source File

SOURCE=..\exinq.c
# End Source File
# Begin Source File

SOURCE=..\exopen.c
# End Source File
# Begin Source File

SOURCE=..\exopts.c
# End Source File
# Begin Source File

SOURCE=..\expclb.c
# End Source File
# Begin Source File

SOURCE=..\expcns.c
# End Source File
# Begin Source File

SOURCE=..\expcon.c
# End Source File
# Begin Source File

SOURCE=..\expcor.c
# End Source File
# Begin Source File

SOURCE=..\expcss.c
# End Source File
# Begin Source File

SOURCE=..\expeat.c
# End Source File
# Begin Source File

SOURCE=..\expelb.c
# End Source File
# Begin Source File

SOURCE=..\expelc.c
# End Source File
# Begin Source File

SOURCE=..\expem.c
# End Source File
# Begin Source File

SOURCE=..\expenm.c
# End Source File
# Begin Source File

SOURCE=..\expev.c
# End Source File
# Begin Source File

SOURCE=..\expfrm.c
# End Source File
# Begin Source File

SOURCE=..\expgv.c
# End Source File
# Begin Source File

SOURCE=..\expinf.c
# End Source File
# Begin Source File

SOURCE=..\expini.c
# End Source File
# Begin Source File

SOURCE=..\expmap.c
# End Source File
# Begin Source File

SOURCE=..\expmp.c
# End Source File
# Begin Source File

SOURCE=..\expnm.c
# End Source File
# Begin Source File

SOURCE=..\expnnm.c
# End Source File
# Begin Source File

SOURCE=..\expnp.c
# End Source File
# Begin Source File

SOURCE=..\expns.c
# End Source File
# Begin Source File

SOURCE=..\expnsd.c
# End Source File
# Begin Source File

SOURCE=..\expnv.c
# End Source File
# Begin Source File

SOURCE=..\expp.c
# End Source File
# Begin Source File

SOURCE=..\exppa.c
# End Source File
# Begin Source File

SOURCE=..\exppn.c
# End Source File
# Begin Source File

SOURCE=..\expqa.c
# End Source File
# Begin Source File

SOURCE=..\expsp.c
# End Source File
# Begin Source File

SOURCE=..\expss.c
# End Source File
# Begin Source File

SOURCE=..\expssd.c
# End Source File
# Begin Source File

SOURCE=..\exptim.c
# End Source File
# Begin Source File

SOURCE=..\expvan.c
# End Source File
# Begin Source File

SOURCE=..\expvnm.c
# End Source File
# Begin Source File

SOURCE=..\expvp.c
# End Source File
# Begin Source File

SOURCE=..\expvpc.c
# End Source File
# Begin Source File

SOURCE=..\expvtt.c
# End Source File
# Begin Source File

SOURCE=..\exupda.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\exodusII.h
# End Source File
# Begin Source File

SOURCE=..\..\include\exodusII_int.h
# End Source File
# End Group
# End Target
# End Project
