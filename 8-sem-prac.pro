TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main/main.c \
    src/common/vectors.c \
    src/kernel/command_line_parser.c \
    src/kernel/general_solver_data.c \
    src/sparse/msr_matrix.c \
    src/common/math_utils.c \
    src/linear_ops/vector_ops.c \
    src/kernel/cgs_solver.c

HEADERS += \
    src/common/debug_utils.h \
    src/common/vectors.h \
    src/common/vectors_fwd.h \
    src/kernel/command_line_parser.h \
    src/kernel/kernel_typedefs.h \
    src/kernel/general_solver_data.h \
    src/sparse/msr_matrix.h \
    src/common/math_utils.h \
    src/linear_ops/vector_ops.h \
    src/kernel/cgs_solver.h \
    src/sparse/sparse_base_format.h

INCLUDEPATH += src

QMAKE_CFLAGS = -std=c89 -ansi -Wall -Werror -Wextra -pedantic-errors
QMAKE_CFLAGS_DEBUG += -DDEBUG

LIBS += -lm -lpthread
