TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main/main.c \
    src/common/vectors.c \
    src/kernel/command_line_parser.c \
    src/sparse/msr_matrix.c \
    src/common/math_utils.c \
    src/linear_ops/vector_ops.c \
    src/kernel/cgs_solver.c \
    src/kernel/input/rhs.c \
    src/kernel/central_diff_solver.c \
    src/kernel/solver_core_workspace.c \
    src/kernel/input/t0_functions.c \
    src/sparse/sparse_base_format.c \
    src/kernel/input/test_solutions.c

HEADERS += \
    src/common/debug_utils.h \
    src/common/vectors.h \
    src/common/vectors_fwd.h \
    src/kernel/command_line_parser.h \
    src/kernel/kernel_typedefs.h \
    src/sparse/msr_matrix.h \
    src/common/math_utils.h \
    src/linear_ops/vector_ops.h \
    src/kernel/cgs_solver.h \
    src/sparse/sparse_base_format.h \
    src/kernel/cgs_solver_private.h \
    src/kernel/command_line_parser_private.h \
    src/kernel/central_diff_solver.h \
    src/kernel/input/rhs.h \
    src/kernel/central_diff_solver_private.h \
    src/kernel/solver_core_workspace.h \
    src/kernel/solver_core_workspace_private.h \
    src/kernel/input/t0_functions.h \
    src/kernel/input/test_solutions.h

INCLUDEPATH += src

QMAKE_CFLAGS = -std=c89 -ansi -Wall -Werror -Wextra -pedantic-errors
QMAKE_CFLAGS_DEBUG += -DDEBUG

LIBS += -lm -lpthread
