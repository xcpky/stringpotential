//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
const std = @import("std");
const testing = std.testing;
const log = std.math.log;
const pow = std.math.pow;
const pi = std.math.pi;
const gsl = @cImport({
    //@cDefine("HAVE_INLINE", .{});
    @cInclude("gsl/gsl_linalg.h");
    @cInclude("gsl/gsl_complex.h");
    @cInclude("gsl/gsl_complex_math.h");
    @cInclude("gsl/gsl_matrix_complex_double.h");
    @cInclude("gsl/gsl_integration.h");
});
const complex = gsl.gsl_complex;
const matrix = gsl.gsl_matrix_complex;
const inverse = gsl.gsl_complex_inverse;
const add = gsl.gsl_complex_add;
const mul = gsl.gsl_complex_mul;
const mul_real = gsl.gsl_complex_mul_real;
const sub = gsl.gsl_complex_sub;
const comp = gsl.gsl_complex_rect;

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

const App = struct {
    Ngauss: usize,
    Lambda: f64,
    epsilon: f64,
    T: *matrix,
    V: *matrix,
    G: *matrix,
    xi: []f64,
    wi: []f64,
    allocator: std.mem.Allocator,
    table: *gsl.gsl_integration_glfixed_table,
    fn init(allocator: std.mem.Allocator, Ngauss: usize, Lambda: f64, epsilon: f64) !void {
        const n = 2 * (Ngauss + 1);
        const Tmat = gsl.gsl_matrix_complex_alloc(n, n) orelse return error.OutOfMemory;
        const Vmat = gsl.gsl_matrix_complex_alloc(n, n) orelse return error.OutOfMemory;
        const Gmat = gsl.gsl_matrix_complex_alloc(n, n) orelse return error.OutOfMemory;
        const table = gsl.gsl_integration_glfixed_table_alloc(Ngauss);
        const xi: []f64 = allocator.alloc(f64, Ngauss) orelse return error.OutOfMemory;
        const wi: []f64 = allocator.alloc(f64, Ngauss) orelse return error.OutOfMemory;
        for (0..Ngauss) |i| {
            gsl.gsl_integration_glfixed_point(0, Lambda, &xi[i], &wi[i], table);
        }
        return App{
            Ngauss,
            Lambda,
            epsilon,
            Tmat,
            Vmat,
            Gmat,
            xi,
            wi,
            allocator,
            table,
        };
    }

    fn deinit(self: *App) !void {
        gsl.gsl_matrix_complex_free(self.T);
        gsl.gsl_matrix_complex_free(self.V);
        gsl.gsl_matrix_complex_free(self.G);
        gsl.gsl_integration_glfixed_table_free(self.table);
        self.allocator.free(self.xi);
        self.allocator.free(self.wi);
    }

    fn gmat(self: *App, E: f64) !void {
        inline for (0..2) |i| {
            const dE = E - delta[i];
            const data = self.G.data;
            const mU = mu[i];
            const x0 = std.math.sqrt(2 * mU * dE);
            var int = comp(0, 0);
            for (self.xi, self.wi) |x, w| {
                var new = comp(E - square(x) / 2 / mU, self.epsilon);
                new = inverse(new);
                new = mul_real(new, w);
                int = add(int, new);
            }
            var tmp = comp(mU * x0 * (log((self.Lambda + x0) / (self.Lambda - x0)) - pi), 0);
            int = mul_real(int, square(x0));
            tmp = sub(tmp, int);
            const tda = self.G.tda;
            data[self.Ngauss * self.Ngauss] = mul_real(tmp, 1 / square(pi) / 2);
            for (i * (self.Ngauss + 1)..i * self.Ngauss + self.Ngauss) |m| {
                data[m * tda + m] = mul_real(inverse(comp(dE - square(self.xi[m]) / 2 / mU, self.epsilon)), square(self.xi[m]) * self.wi[m] / 2 / square(pi));
            }
        }
    }

    fn vmat(self: *App, E: f64) !void {
        const x = std.ArrayList(f64).fromOwnedSlice(self.allocator, self.xi);
        x.addOne();
        inline for (0..2) |i| {
            const dE = E - delta[i];
            const data = self.V.data;
            const mU = mu[i];
            const x0 = std.math.sqrt(2*mU*dE);
            x[-1] = x0;
        }
    }
};

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

fn Omega(comptime alpha: u8, comptime beta: u8) type {
    return struct {
        pub fn compute(p: f64, pprime: f64) f64 {
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
    };
}

fn Omega_prime(comptime alpha: u8, comptime beta: u8) type {
    return struct {
        pub fn compute(p: f64, pprime: f64) f64 {
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
    };
}

fn O(comptime alpha: u8, comptime beta: u8, E: f64, p: f64, pprime: f64, m: f64) f64 {
    const omega = Omega(alpha, beta).compute(p, pprime);
    const omega_prime = Omega_prime(alpha, beta).compute(p, pprime);
    const tmp = []f64{ E - (m + square(p - pprime) / 2 / m), E - (m + square(p + pprime) / 2 / m) };
    return -1 / p / pprime / 4 * (log(f64, (tmp[0] - omega) / (tmp[1] - omega)) + log(f64, (tmp[0] - omega_prime) / (tmp[1] - omega_prime)));
}

fn V_OME(comptime alpha: u8, comptime beta: u8) type {
    return struct {
        pub fn compute(E: f64, p: f64, pprime: f64) f64 {
            return switch (alpha) {
                0 => switch (beta) {
                    0 => -3 * (3 * O(0, 0, E, p, pprime, m_pi) + O(0, 0, E, p, pprime, m_eta) / 3),
                    1 => pow(f64, 2, 3 / 2) * O(0, 1, E, p, pprime, m_K),
                    else => @compileError("Invalid beta"),
                },
                1 => switch (beta) {
                    0 => pow(f64, 2, 3 / 2) * O(1, 0, E, pprime, m_K),
                    1 => 2 / 3 * O(1, 1, E, p, pprime, m_eta),
                    else => @compileError("Invalid beta"),
                },
                else => @compileError("Invalid alpha"),
            };
        }
    };
}

test "basic add functionality" {
    try testing.expect(add(3, 7) == 10);
}
