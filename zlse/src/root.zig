//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
const std = @import("std");
const testing = std.testing;
const log = std.math.log;
const pow = std.math.pow;
const pi = std.math.pi;
const debug = false;

const gsl = @cImport({
    // @cDefine("HAVE_INLINE", {});
    @cInclude("gsl/gsl_linalg.h");
    @cInclude("gsl/gsl_complex.h");
    // @cInclude("gsl/gsl_permute_vector_complex_long_double.h");
    @cInclude("gsl/gsl_complex_math.h");
    @cInclude("gsl/gsl_matrix_complex_double.h");
    @cInclude("gsl/gsl_integration.h");
    @cInclude("gsl/gsl_blas.h");
    @cInclude("gsl/gsl_cblas.h");
});
const complex = gsl.gsl_complex;
const matrix = gsl.gsl_matrix_complex;
const inverse = gsl.gsl_complex_inverse;
const add = gsl.gsl_complex_add;
const mul = gsl.gsl_complex_mul;
const mul_real = gsl.gsl_complex_mul_real;
const sub = gsl.gsl_complex_sub;
const comp = gsl.gsl_complex_rect;
const matrix_alloc = gsl.gsl_matrix_complex_alloc;
const matrix_free = gsl.gsl_matrix_complex_free;
const matrix_set_zero = gsl.gsl_matrix_complex_set_zero;
const matrix_set = gsl.gsl_matrix_complex_set;
const matrix_memcpy = gsl.gsl_matrix_complex_memcpy;
const matrix_scale = gsl.gsl_matrix_complex_scale;
const matrix_get = gsl.gsl_matrix_complex_get;
const CBLAS_ORDER = gsl.CBLAS_ORDER;
const CblasColMajor = gsl.CblasColMajor;
const CblasRowMajor = gsl.CblasRowMajor;
const CblasNoTrans = gsl.CblasNoTrans;

const f_pi: f64 = 0.092;
const g_pi: f64 = 0.5704;
const m_pi: f64 = 0.138039407;
const m_K: f64 = 0.498;
const m_eta: f64 = 0.548;
const m_eta_s: f64 = 0.690625803;
const m_B: f64 = 5.27934;
const m_B_star: f64 = 5.32471;
const m_B_s: f64 = 5.36692;
const m_B_star_s: f64 = 5.4154;
const m11: f64 = m_B;
const m12: f64 = m_B_star;
const m21: f64 = m_B_s;
const m22: f64 = m_B_star_s;
const delta = [_]f64{ 0, m21 + m22 - m11 - m12 };
const mu = [_]f64{ m11 * m12 / (m11 + m12), m21 * m22 / (m21 + m22) };

inline fn square(x: f64) f64 {
    return x * x;
}

// inline fn comp(x: f64, y: f64) complex {
//     return complex{
//         .dat = [2]f64{x, y},
//     };
// }
//
const LSE = struct {
    Ngauss: usize,
    Lambda: f64,
    epsilon: f64,
    E: f64,
    T: *matrix,
    V: *matrix,
    G: *matrix,
    xi: []f64,
    wi: []f64,
    x0: [2]f64,
    allocator: std.mem.Allocator,
    table: *gsl.gsl_integration_glfixed_table,
    fn init(allocator: std.mem.Allocator, Ngauss: usize, Lambda: f64, epsilon: f64, E: f64) !LSE {
        const n = 2 * (Ngauss + 1);
        const Tmat = matrix_alloc(n, n) orelse return error.OutOfMemory;
        const Vmat = matrix_alloc(n, n) orelse return error.OutOfMemory;
        const Gmat = matrix_alloc(n, n) orelse return error.OutOfMemory;
        // if (debug == true) {
        //     std.debug.print("in initialization\nn: {}\n", .{n});
        //     std.debug.print("gmat: size1 {}, size2 {}\n", .{ Gmat.*.size1, Gmat.*.size2 });
        //     std.debug.print("vmat: size1 {}, size2 {}\n", .{ Vmat.*.size1, Vmat.*.size2 });
        //     std.debug.print("tmat: size1 {}, size2 {}\n", .{ Tmat.*.size1, Tmat.*.size2 });
        // }
        const table = gsl.gsl_integration_glfixed_table_alloc(Ngauss);
        const xi: []f64 = try allocator.alloc(f64, Ngauss);
        const wi: []f64 = try allocator.alloc(f64, Ngauss);
        for (0..Ngauss) |i| {
            _ = gsl.gsl_integration_glfixed_point(0, Lambda, i, &xi[i], &wi[i], table);
        }
        var x0: [2]f64 = undefined;
        inline for (0..2) |i| {
            const dE = E - delta[i];
            const mU = mu[i];
            const tmp = std.math.sqrt(2 * mU * dE);
            x0[i] = tmp;
        }
        return LSE{
            .Ngauss = Ngauss,
            .Lambda = Lambda,
            .epsilon = epsilon,
            .E = E,
            .T = Tmat,
            .V = Vmat,
            .G = Gmat,
            .xi = xi,
            .wi = wi,
            .x0 = [2]f64{ x0[0], x0[1] },
            .allocator = allocator,
            .table = table,
        };
    }

    fn refresh(self: *LSE, Lambda: f64, E: f64) void {
        self.Lambda = Lambda;
        self.E = E;
        inline for (0..2) |i| {
            const dE = E - delta[i];
            const mU = mu[i];
            const tmp = std.math.sqrt(2 * mU * dE);
            self.x0[i] = tmp;
        }
    }

    fn deinit(self: *LSE) void {
        matrix_free(self.T);
        matrix_free(self.V);
        matrix_free(self.G);
        gsl.gsl_integration_glfixed_table_free(self.table);
        self.allocator.free(self.xi);
        self.allocator.free(self.wi);
    }

    fn gmat(self: *LSE) !void {
        matrix_set_zero(self.G);
        inline for (0..2) |i| {
            const dE = self.E - delta[i];
            const mU = mu[i];
            const x0 = self.x0[i];
            var int = comp(0, 0);
            std.debug.print("dE: {e:.3} mU: {e:.3} epsilon: {e:.3}\n", .{ dE, mU, self.epsilon });
            for (self.xi, self.wi) |x, w| {
                var new = comp(self.E - square(x) / 2 / mU, self.epsilon);
                new = inverse(new);
                new = mul_real(new, w);
                int = add(int, new);
            }
            var tmp = comp(mU * x0 * (log(f64, std.math.e, (self.Lambda + x0) / (self.Lambda - x0)) - pi), 0);
            int = mul_real(int, square(x0));
            tmp = sub(tmp, int);
            // const tda = self.G.tda;
            // data[self.Ngauss * self.Ngauss] = mul_real(tmp, 1 / square(pi) / 2);
            const ii = self.Ngauss + i * (self.Ngauss + 1);
            matrix_set(self.G, ii, ii, mul_real(tmp, 1 / square(pi) / 2));
            // std.debug.print("ii: {} {d:.3} {d:.3}\n", .{ii, self.G.data[2*(ii*tda + ii)], self.G.data[2*(ii*tda + ii) + 1]});
            for (0..self.Ngauss) |m| {
                const pos = m + i * (self.Ngauss + 1);
                const ele = mul_real(inverse(comp(dE - square(self.xi[m]) / 2 / mU, self.epsilon)), square(self.xi[m]) * self.wi[m] / 2 / square(pi));
                const lhs = inverse(comp(dE - square(self.xi[m]) / 2 / mU, self.epsilon));
                const rhs = square(self.xi[m]) * self.wi[m] / 2 / square(pi);
                const denominator = comp(dE - square(self.xi[m]) / 2 / mU, self.epsilon);
                const inv = inverse(denominator);
                std.debug.print("xi: {e:.3}\n", .{self.xi[m]});
                std.debug.print("denominator: {e:.3} + im{e:.3}        inverse: {e:.3} + im{e:.3}\n", .{ denominator.dat[0], denominator.dat[1], inv.dat[0], inv.dat[1] });
                std.debug.print("lhs: {e:.3} + im{e:.3}, rhs: {e:.3}\n\n\n", .{ lhs.dat[0], lhs.dat[1], rhs });
                matrix_set(self.G, pos, pos, ele);
                // std.debug.print("m: {} {e} {e}\n", .{m, self.xi[m], mul_real(inverse(comp(dE - square(self.xi[m]) / 2 / mU, self.epsilon)), square(self.xi[m]) * self.wi[m] / 2 / square(pi)).dat[0]});
            }
        }
    }

    fn vmat(self: *LSE) !void {
        matrix_set_zero(self.V);
        const x: [][]f64 = try self.allocator.alloc([]f64, 2);
        defer self.allocator.free(x);
        for (x) |*arr| {
            arr.* = try self.allocator.alloc(f64, self.Ngauss + 1);
            std.mem.copyForwards(f64, arr.*, self.xi);
        }
        defer for (x) |arr| {
            self.allocator.free(arr);
        };
        inline for (0..2) |i| {
            x[i][self.Ngauss] = self.x0[i];
        }
        // const tda = self.V.tda;
        // const data = self.V.data;
        inline for (0..2) |alpha| {
            inline for (0..2) |beta| {
                for (x[alpha], 0..) |p, idx| {
                    for (x[beta], 0..) |pprime, jdx| {
                        const i = idx + alpha * (self.Ngauss + 1);
                        const j = jdx + beta * (self.Ngauss + 1);
                        // data[2 * (i * tda + j)] = V_OME(alpha, beta, self.E, p, pprime);
                        matrix_set(self.V, i, j, comp(V_OME(alpha, beta, self.E, p, pprime), 0));
                    }
                }
            }
        }
    }

    fn tmat(self: *LSE) !void {
        const n = 2 * (self.Ngauss + 1);

        // Step 1: Compute VG = V * G
        const VG: *matrix = matrix_alloc(n, n) orelse return error.OutOfMemory;
        defer matrix_free(VG);
        matrix_set_zero(VG);

        // Perform matrix multiplication: VG = V * G
        const alpha = comp(1.0, 0.0);
        // if (debug == true) {
        //     std.debug.print("in tmat\n", .{});
        //     std.debug.print("gmat: size1 {}, size2 {}, tda{}\n", .{ self.G.size1, self.G.size2, self.G.tda });
        //     std.debug.print("vmat: size1 {}, size2 {}, tda{}\n", .{ self.V.size1, self.V.size2, self.V.tda });
        //     std.debug.print("tmat: size1 {}, size2 {}, tda{}\n", .{ self.T.size1, self.T.size2, self.T.tda });
        //     std.debug.print("vg: size1 {}, size2 {}, tda{}\n", .{ VG.size1, VG.size2, VG.tda });
        // }
        const beta = comp(0.0, 0.0);
        gsl.cblas_zgemm(
            CblasRowMajor, // Matrix storage order
            CblasNoTrans, // No transpose A
            CblasNoTrans, // No transpose B
            @intCast(n), // Rows of A (mat1) and C (mat3)
            @intCast(n), // Columns of B (mat2) and C (mat3)
            @intCast(n), // Columns of A / rows of B (inner dimension)
            &alpha.dat[0], // α as [real, imag] pointer
            self.V.data, // A data (column-major)
            @intCast(self.V.tda), // Leading dimension of A
            self.G.data, // B data (column-major)
            @intCast(self.G.tda), // Leading dimension of B
            &beta.dat[0], // β as [real, imag] pointer
            VG.data, // C data (column-major)
            @intCast(VG.tda), // Leading dimension of C
        );
        // if (gsl.gsl_blas_zgemm(
        //     gsl.CblasNoTrans,
        //     gsl.CblasNoTrans,
        //     alpha,
        //     self.V,
        //     self.G,
        //     beta,
        //     VG,
        // ) != gsl.GSL_SUCCESS) {
        //     return error.MatrixMultiplyFailed;
        // }
        //
        // Step 2: Compute I - VG
        const I_minus_VG = matrix_alloc(n, n) orelse return error.OutOfMemory;
        defer matrix_free(I_minus_VG);
        _ = matrix_memcpy(I_minus_VG, VG);
        _ = matrix_scale(I_minus_VG, comp(-1.0, 0.0)); // I_minus_VG = -VG

        // Add identity matrix: I_minus_VG = I - VG
        for (0..n) |i| {
            const diag = matrix_get(I_minus_VG, i, i);
            const one = comp(1.0, 0.0);
            matrix_set(I_minus_VG, i, i, add(diag, one));
        }

        // Step 3: Invert (I - VG) using LU decomposition
        const perm = gsl.gsl_permutation_alloc(n) orelse return error.OutOfMemory;
        defer gsl.gsl_permutation_free(perm);

        var signum: c_int = undefined;
        if (gsl.gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) != gsl.GSL_SUCCESS) {
            return error.LUDecompFailed;
        }

        const inv_I_minus_VG: *matrix = matrix_alloc(n, n) orelse return error.OutOfMemory;
        defer matrix_free(inv_I_minus_VG);
        if (gsl.gsl_linalg_complex_LU_invert(I_minus_VG, perm, inv_I_minus_VG) != gsl.GSL_SUCCESS) {
            return error.MatrixInverseFailed;
        }

        gsl.cblas_zgemm(
            CblasRowMajor, // Matrix storage order
            CblasNoTrans, // No transpose A
            CblasNoTrans, // No transpose B
            @intCast(n), // Rows of A (mat1) and C (mat3)
            @intCast(n), // Columns of B (mat2) and C (mat3)
            @intCast(n), // Columns of A / rows of B (inner dimension)
            &alpha.dat[0], // α as [real, imag] pointer
            inv_I_minus_VG.data, // A data (column-major)
            @intCast(inv_I_minus_VG.tda), // Leading dimension of A
            self.V.data, // B data (column-major)
            @intCast(self.V.tda), // Leading dimension of B
            &beta.dat[0], // β as [real, imag] pointer
            self.T.data, // C data (column-major)
            @intCast(self.T.tda), // Leading dimension of C
        );
        // Step 4: Compute T = inv(I - VG) * V
        // if (gsl.gsl_blas_zgemm(
        //     gsl.CblasNoTrans,
        //     gsl.CblasNoTrans,
        //     alpha,
        //     inv_I_minus_VG,
        //     self.V,
        //     beta,
        //     self.T,
        // ) != gsl.GSL_SUCCESS) {
        //     return error.MatrixMultiplyFailed;
        // }
    }

    fn run(self: *LSE) !void {
        try self.gmat();
        try self.vmat();
        try self.tmat();
    }

    fn compute(self: *LSE, Lambda: f64, E: f64) !void {
        self.refresh(Lambda, E);
        try self.run();
    }
};

// Export the App type as an opaque handle
pub const LSEHandle = *LSE;

// Create a new App instance
pub export fn lse_new(
    Ngauss: u64,
    Lambda: f64,
    epsilon: f64,
    E: f64,
) ?*LSE {
    const allocator = std.heap.c_allocator;
    const app = allocator.create(LSE) catch return null;
    app.* = LSE.init(allocator, Ngauss, Lambda, epsilon, E) catch {
        allocator.destroy(app);
        return null;
    };
    return app;
}

// Run the LSE solver
pub export fn lse_run(app: *LSE) i32 {
    app.run() catch return -1;
    return 0;
}

pub export fn lse_compute(app: *LSE, Lambda: f64, E: f64) void {
    app.compute(Lambda, E) catch @panic("shit");
}

// Free the App instance
pub export fn lse_deinit(app: *LSE) void {
    app.deinit();
    const allocator = app.allocator;
    allocator.destroy(app);
}

// Get pointers to matrix data (G, V, T)
pub export fn lse_get_g_data(app: *LSE) [*]f64 {
    return @ptrCast(@alignCast(app.G.data));
}

pub export fn lse_get_g_size(app: *LSE, rows: *c_uint, cols: *c_uint) void {
    rows.* = @intCast(app.G.size1);
    cols.* = @intCast(app.G.size2);
}

pub export fn lse_get_v_data(app: *LSE) [*]f64 {
    return @ptrCast(@alignCast(app.V.data));
}

pub export fn lse_get_v_size(app: *LSE, rows: *c_uint, cols: *c_uint) void {
    rows.* = @intCast(app.V.size1);
    cols.* = @intCast(app.V.size2);
}

pub export fn lse_get_t_data(app: *LSE) [*]f64 {
    return @ptrCast(@alignCast(app.T.data));
}

pub export fn lse_get_t_size(app: *LSE, rows: *c_uint, cols: *c_uint) void {
    rows.* = @intCast(app.T.size1);
    cols.* = @intCast(app.T.size2);
}

fn integrate(
    comptime Context: type,
    ctx: Context,
    comptime f: fn (Context, f64) f64,
    nodes: []const f64,
    weights: []const f64,
) complex {
    var sum: complex = comp(0, 0);
    for (nodes, weights) |x, w| {
        sum = add(f(ctx, x) * w, sum);
    }
    return sum;
}

fn Omega(comptime alpha: u8, comptime beta: u8, p: f64, pprime: f64) f64 {
    return switch (alpha) {
        0 => switch (beta) {
            0 => 2 * m_B + (p * p + pprime * pprime) / 2 / m_B,
            1 => m_B + pprime * pprime / 2 / m_B + m_B_s + p * p / 2 / m_B_s,
            else => @compileError("Invalid beta"),
        },
        1 => switch (beta) {
            0 => m_B_s + pprime * pprime / 2 / m_B_s + m_B + p * p / 2 / m_B,
            1 => 2 * m_B_s + (p * p + pprime * pprime) / 2 / m_B_s,
            else => @compileError("Invalid beta"),
        },
        else => @compileError("Invalid alpha"),
    };
}

fn Omega_prime(comptime alpha: u8, comptime beta: u8, p: f64, pprime: f64) f64 {
    return switch (alpha) {
        0 => switch (beta) {
            0 => 2 * m_B_star + (p * p + pprime * pprime) / 2 / m_B_star,
            1 => m_B_star + pprime * pprime / 2 / m_B_star + m_B_star_s + p * p / 2 / m_B_star_s,
            else => @compileError("Invalid beta"),
        },
        1 => switch (beta) {
            0 => m_B_star_s + pprime * pprime / 2 / m_B_star_s + m_B_star + p * p / 2 / m_B_star,
            1 => 2 * m_B_star_s + (p * p + pprime * pprime) / 2 / m_B_star_s,
            else => @compileError("Invalid beta"),
        },
        else => @compileError("Invalid alpha"),
    };
}

fn O(comptime alpha: u8, comptime beta: u8, E: f64, p: f64, pprime: f64, m: f64) f64 {
    const omega = Omega(alpha, beta, p, pprime);
    const omega_prime = Omega_prime(alpha, beta, p, pprime);
    const tmp = [2]f64{ E - (m + square(p - pprime) / 2 / m), E - (m + square(p + pprime) / 2 / m) };
    return -1 / p / pprime / 4 * (log(f64, std.math.e, (tmp[0] - omega) / (tmp[1] - omega)) + log(f64, std.math.e, (tmp[0] - omega_prime) / (tmp[1] - omega_prime)));
}

fn V_OME(comptime alpha: u8, comptime beta: u8, E: f64, p: f64, pprime: f64) f64 {
    return switch (alpha) {
        0 => switch (beta) {
            0 => -3 * (3 * O(0, 0, E, p, pprime, m_pi) + O(0, 0, E, p, pprime, m_eta) / 3),
            1 => pow(f64, 2, 3 / 2) * O(0, 1, E, p, pprime, m_K),
            else => @compileError("Invalid beta"),
        },
        1 => switch (beta) {
            0 => pow(f64, 2, 3 / 2) * O(1, 0, E, p, pprime, m_K),
            1 => 2 / 3 * O(1, 1, E, p, pprime, m_eta),
            else => @compileError("Invalid beta"),
        },
        else => @compileError("Invalid alpha"),
    };
}

fn test_blas(comptime n: u8) !void {
    const mat1: *matrix = matrix_alloc(n, n) orelse return error.OutOfMemory;
    const mat2: *matrix = matrix_alloc(n, n) orelse return error.OutOfMemory;
    const mat3: *matrix = matrix_alloc(n, n) orelse return error.OutOfMemory;
    defer matrix_free(mat1);
    defer matrix_free(mat2);
    defer matrix_free(mat3);
    gsl.gsl_matrix_complex_set_zero(mat1);
    gsl.gsl_matrix_complex_set_zero(mat2);
    matrix_set(mat1, 0, 1, comp(1, 1));
    matrix_set(mat1, 1, 0, comp(-1, 0));
    matrix_set(mat2, 0, 0, comp(2, 0));
    matrix_set(mat2, 0, 1, comp(3, -1));
    matrix_set(mat2, 1, 0, comp(1, -1));
    inline for (0..2) |i| {
        inline for (0..2) |j| {
            const pos = i * 2 + j;
            const val = mat1.data[2 * pos];
            // const val = matrix_get(mat2, i, j);
            std.debug.print("{d:.3} + Im{d:.3}  ", .{ val, mat1.data[2 * pos + 1] });
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("\n", .{});
    inline for (0..2) |i| {
        inline for (0..2) |j| {
            const pos = i * 2 + j;
            const val = mat2.data[2 * pos];
            // const val = matrix_get(mat2, i, j);
            std.debug.print("{d:.3} + Im{d:.3}  ", .{ val, mat2.data[2 * pos + 1] });
        }
        std.debug.print("\n", .{});
    }
    matrix_set_zero(mat3);
    const alpha = comp(1, 0);
    std.debug.print("alpha: {d:.3} + Im{d:.3}\n", .{ alpha.dat[0], alpha.dat[1] });
    const beta = comp(0, 0);
    // Verify tda matches row count
    std.debug.assert(mat1.tda == mat1.size1);
    std.debug.assert(mat2.tda == mat2.size1);
    std.debug.assert(mat3.tda == mat3.size1);
    // Verify dimensions
    std.debug.assert(mat1.size1 == n and mat1.size2 == n);
    std.debug.assert(mat2.size1 == n and mat2.size2 == n);
    std.debug.assert(mat3.size1 == n and mat3.size2 == n);
    // _ = gsl.gsl_blas_zgemm(gsl.CblasNoTrans, gsl.CblasNoTrans, alpha, mat1, mat2, beta, mat3);
    gsl.cblas_zgemm(
        CblasRowMajor, // Matrix storage order
        CblasNoTrans, // No transpose A
        CblasNoTrans, // No transpose B
        n, // Rows of A (mat1) and C (mat3)
        n, // Columns of B (mat2) and C (mat3)
        n, // Columns of A / rows of B (inner dimension)
        &alpha.dat[0], // α as [real, imag] pointer
        mat1.data, // A data (column-major)
        @intCast(mat1.tda), // Leading dimension of A
        mat2.data, // B data (column-major)
        @intCast(mat2.tda), // Leading dimension of B
        &beta.dat[0], // β as [real, imag] pointer
        mat3.data, // C data (column-major)
        @intCast(mat3.tda), // Leading dimension of C
    );
    // const data = mat3.data;
    // for (mat2.data) |val| {
    //     std.debug.print("{d:.3} ", .{val});
    // }
    // std.debug.print("\n", .{});
    inline for (0..2) |i| {
        inline for (0..2) |j| {
            // const pos = i*2 + j;
            // const val = mat3.data[2*pos];
            const val = matrix_get(mat3, i, j);
            std.debug.print("{d:.3} + Im{d:.3}  ", .{val.dat[0], val.dat[1]});
        }
        std.debug.print("\n", .{});
    }
}

fn test_mat() !void {
    var x = comp(1, 2);
    x = inverse(x);
    std.debug.print("{e:.3} + Im{e:.3}\n\n\n\n", .{ x.dat[0], x.dat[1] });
}

test "small test" {}

test "basic functionality" {
    // try testing.expect(add(3, 7) == 10);
    // try test_mat();
    // var lse = try LSE.init(std.heap.c_allocator, 2, 4, 0.000001, 0.2);
    // defer lse.deinit();
    // try lse.run();
    // const n = 2*(lse.Ngauss + 1);
    // for (0..n) |i| {
    //     for (0..n) |j| {
    //         // const pos = i*n + j;
    //         const val = matrix_get(lse.G, i, j);
    //         std.debug.print("{d:.5} ", .{val.dat[1]});
    //     }
    //     std.debug.print("\n", .{});
    // }
    try test_blas(2);
    // std.debug.print("gmat: size1 {}, size2 {}\n", .{ lse.G.size1, lse.G.size2 });
    // std.debug.print("vmat: size1 {}, size2 {}\n", .{ lse.V.size1, lse.V.size2 });
    // std.debug.print("tmat: size1 {}, size2 {}\n", .{ lse.T.size1, lse.T.size2 });
    // try lse.run();
}
