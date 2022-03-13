cmake_minimum_required(VERSION 2.8)
project(INSTALL C Fortran)

add_executable(xintface Fintface.f Cintface.c)
set_property(TARGET xintface PROPERTY POSITION_INDEPENDENT_CODE TRUE)
