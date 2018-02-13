TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main/main.c \
    src/common/vectors.c \
    src/kernel/central_differences_solver.c \
    src/kernel/command_line_parser.c

HEADERS += \
    src/common/debug_utils.h \
    src/common/vectors.h \
    src/kernel/central_differences_solver.h \
    src/common/vectors_fwd.h \
    src/kernel/command_line_parser.h \
    src/kernel/kernel_typedefs.h

INCLUDEPATH += src
