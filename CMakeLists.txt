project(fem-cloth)
find_package(Qt4 COMPONENTS QtCore QtGui QtOpenGL QtXml REQUIRED)
find_package(OpenGL REQUIRED)
include_directories(/usr/include/eigen3)
include(${QT_USE_FILE})
cmake_minimum_required(VERSION 2.8)
add_definitions(-std=c++11)
set_source_files_properties(triangle.c PROPERTIES COMPILE_DEFINITIONS "TRILIBRARY;REDUCED")
QT4_WRAP_CPP(SRCS Viewer.hpp Simulator.hpp)
add_executable(${PROJECT_NAME}
main.cpp ${SRCS}
triangle.h triangle.c
Viewer.hpp Viewer.cpp
Simulator.hpp Simulator.cpp
)

target_link_libraries(${PROJECT_NAME} ${QT_LIBRARIES} qglviewer-qt4 GLEW
${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

