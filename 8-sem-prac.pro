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
    src/kernel/input/test_solutions.c \
    src/3rd_party/laspack/eigenval.c \
    src/3rd_party/laspack/errhandl.c \
    src/3rd_party/laspack/factor.c \
    src/3rd_party/laspack/itersolv.c \
    src/3rd_party/laspack/matrix.c \
    src/3rd_party/laspack/mlsolv.c \
    src/3rd_party/laspack/operats.c \
    src/3rd_party/laspack/precond.c \
    src/3rd_party/laspack/qmatrix.c \
    src/3rd_party/laspack/rtc.c \
    src/3rd_party/laspack/vector.c \
    src/sparse/laspack_matrix.c \
    src/sparse/laspack_vector.c \
    src/io/table_io.c \
    src/kernel/solver_tester.c \
    src/io/gnuplot_io.c \
    src/kernel/half_nodes_values.c \
    src/kernel/nodes_values.c \
    src/kernel/linear_system_composer.c \
    src/kernel/sokolov_solver.c \
    src/kernel/mesh_info.c

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
    src/kernel/input/test_solutions.h \
    src/3rd_party/laspack/copyrght.h \
    src/3rd_party/laspack/eigenval.h \
    src/3rd_party/laspack/elcmp.h \
    src/3rd_party/laspack/errhandl.h \
    src/3rd_party/laspack/factor.h \
    src/3rd_party/laspack/itersolv.h \
    src/3rd_party/laspack/lastypes.h \
    src/3rd_party/laspack/matrix.h \
    src/3rd_party/laspack/mlsolv.h \
    src/3rd_party/laspack/operats.h \
    src/3rd_party/laspack/precond.h \
    src/3rd_party/laspack/qmatrix.h \
    src/3rd_party/laspack/rtc.h \
    src/3rd_party/laspack/vector.h \
    src/3rd_party/laspack/version.h \
    src/sparse/laspack_matrix.h \
    src/sparse/laspack_vector.h \
    src/io/table_io.h \
    src/kernel/solver_tester.h \
    src/io/table_io_private.h \
    src/io/gnuplot_io.h \
    src/kernel/half_nodes_values.h \
    src/kernel/nodes_values.h \
    src/kernel/mesh_info.h \
    src/kernel/linear_system_composer.h \
    src/kernel/sokolov_solver.h \
    src/kernel/sokolov_solver_private.h

INCLUDEPATH += src

QMAKE_CFLAGS = -std=c89 -ansi -Wall -Werror -Wextra
QMAKE_CFLAGS_DEBUG += -DDEBUG

LIBS += -lm -lpthread
