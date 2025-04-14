pub const __builtin_bswap16 = @import("std").zig.c_builtins.__builtin_bswap16;
pub const __builtin_bswap32 = @import("std").zig.c_builtins.__builtin_bswap32;
pub const __builtin_bswap64 = @import("std").zig.c_builtins.__builtin_bswap64;
pub const __builtin_signbit = @import("std").zig.c_builtins.__builtin_signbit;
pub const __builtin_signbitf = @import("std").zig.c_builtins.__builtin_signbitf;
pub const __builtin_popcount = @import("std").zig.c_builtins.__builtin_popcount;
pub const __builtin_ctz = @import("std").zig.c_builtins.__builtin_ctz;
pub const __builtin_clz = @import("std").zig.c_builtins.__builtin_clz;
pub const __builtin_sqrt = @import("std").zig.c_builtins.__builtin_sqrt;
pub const __builtin_sqrtf = @import("std").zig.c_builtins.__builtin_sqrtf;
pub const __builtin_sin = @import("std").zig.c_builtins.__builtin_sin;
pub const __builtin_sinf = @import("std").zig.c_builtins.__builtin_sinf;
pub const __builtin_cos = @import("std").zig.c_builtins.__builtin_cos;
pub const __builtin_cosf = @import("std").zig.c_builtins.__builtin_cosf;
pub const __builtin_exp = @import("std").zig.c_builtins.__builtin_exp;
pub const __builtin_expf = @import("std").zig.c_builtins.__builtin_expf;
pub const __builtin_exp2 = @import("std").zig.c_builtins.__builtin_exp2;
pub const __builtin_exp2f = @import("std").zig.c_builtins.__builtin_exp2f;
pub const __builtin_log = @import("std").zig.c_builtins.__builtin_log;
pub const __builtin_logf = @import("std").zig.c_builtins.__builtin_logf;
pub const __builtin_log2 = @import("std").zig.c_builtins.__builtin_log2;
pub const __builtin_log2f = @import("std").zig.c_builtins.__builtin_log2f;
pub const __builtin_log10 = @import("std").zig.c_builtins.__builtin_log10;
pub const __builtin_log10f = @import("std").zig.c_builtins.__builtin_log10f;
pub const __builtin_abs = @import("std").zig.c_builtins.__builtin_abs;
pub const __builtin_labs = @import("std").zig.c_builtins.__builtin_labs;
pub const __builtin_llabs = @import("std").zig.c_builtins.__builtin_llabs;
pub const __builtin_fabs = @import("std").zig.c_builtins.__builtin_fabs;
pub const __builtin_fabsf = @import("std").zig.c_builtins.__builtin_fabsf;
pub const __builtin_floor = @import("std").zig.c_builtins.__builtin_floor;
pub const __builtin_floorf = @import("std").zig.c_builtins.__builtin_floorf;
pub const __builtin_ceil = @import("std").zig.c_builtins.__builtin_ceil;
pub const __builtin_ceilf = @import("std").zig.c_builtins.__builtin_ceilf;
pub const __builtin_trunc = @import("std").zig.c_builtins.__builtin_trunc;
pub const __builtin_truncf = @import("std").zig.c_builtins.__builtin_truncf;
pub const __builtin_round = @import("std").zig.c_builtins.__builtin_round;
pub const __builtin_roundf = @import("std").zig.c_builtins.__builtin_roundf;
pub const __builtin_strlen = @import("std").zig.c_builtins.__builtin_strlen;
pub const __builtin_strcmp = @import("std").zig.c_builtins.__builtin_strcmp;
pub const __builtin_object_size = @import("std").zig.c_builtins.__builtin_object_size;
pub const __builtin___memset_chk = @import("std").zig.c_builtins.__builtin___memset_chk;
pub const __builtin_memset = @import("std").zig.c_builtins.__builtin_memset;
pub const __builtin___memcpy_chk = @import("std").zig.c_builtins.__builtin___memcpy_chk;
pub const __builtin_memcpy = @import("std").zig.c_builtins.__builtin_memcpy;
pub const __builtin_expect = @import("std").zig.c_builtins.__builtin_expect;
pub const __builtin_nanf = @import("std").zig.c_builtins.__builtin_nanf;
pub const __builtin_huge_valf = @import("std").zig.c_builtins.__builtin_huge_valf;
pub const __builtin_inff = @import("std").zig.c_builtins.__builtin_inff;
pub const __builtin_isnan = @import("std").zig.c_builtins.__builtin_isnan;
pub const __builtin_isinf = @import("std").zig.c_builtins.__builtin_isinf;
pub const __builtin_isinf_sign = @import("std").zig.c_builtins.__builtin_isinf_sign;
pub const __has_builtin = @import("std").zig.c_builtins.__has_builtin;
pub const __builtin_assume = @import("std").zig.c_builtins.__builtin_assume;
pub const __builtin_unreachable = @import("std").zig.c_builtins.__builtin_unreachable;
pub const __builtin_constant_p = @import("std").zig.c_builtins.__builtin_constant_p;
pub const __builtin_mul_overflow = @import("std").zig.c_builtins.__builtin_mul_overflow;
pub const wchar_t = c_int;
// /usr/include/bits/floatn.h:83:24: warning: unsupported type: 'Complex'
pub const __cfloat128 = @compileError("unable to resolve typedef child type");
// /usr/include/bits/floatn.h:83:24
pub const _Float128 = f128;
pub const _Float32 = f32;
pub const _Float64 = f64;
pub const _Float32x = f64;
pub const _Float64x = c_longdouble;
pub const div_t = extern struct {
    quot: c_int = @import("std").mem.zeroes(c_int),
    rem: c_int = @import("std").mem.zeroes(c_int),
};
pub const ldiv_t = extern struct {
    quot: c_long = @import("std").mem.zeroes(c_long),
    rem: c_long = @import("std").mem.zeroes(c_long),
};
pub const lldiv_t = extern struct {
    quot: c_longlong = @import("std").mem.zeroes(c_longlong),
    rem: c_longlong = @import("std").mem.zeroes(c_longlong),
};
pub extern fn __ctype_get_mb_cur_max() usize;
pub extern fn atof(__nptr: [*c]const u8) f64;
pub extern fn atoi(__nptr: [*c]const u8) c_int;
pub extern fn atol(__nptr: [*c]const u8) c_long;
pub extern fn atoll(__nptr: [*c]const u8) c_longlong;
pub extern fn strtod(__nptr: [*c]const u8, __endptr: [*c][*c]u8) f64;
pub extern fn strtof(__nptr: [*c]const u8, __endptr: [*c][*c]u8) f32;
pub extern fn strtold(__nptr: [*c]const u8, __endptr: [*c][*c]u8) c_longdouble;
pub extern fn strtol(__nptr: [*c]const u8, __endptr: [*c][*c]u8, __base: c_int) c_long;
pub extern fn strtoul(__nptr: [*c]const u8, __endptr: [*c][*c]u8, __base: c_int) c_ulong;
pub extern fn strtoq(noalias __nptr: [*c]const u8, noalias __endptr: [*c][*c]u8, __base: c_int) c_longlong;
pub extern fn strtouq(noalias __nptr: [*c]const u8, noalias __endptr: [*c][*c]u8, __base: c_int) c_ulonglong;
pub extern fn strtoll(__nptr: [*c]const u8, __endptr: [*c][*c]u8, __base: c_int) c_longlong;
pub extern fn strtoull(__nptr: [*c]const u8, __endptr: [*c][*c]u8, __base: c_int) c_ulonglong;
pub extern fn l64a(__n: c_long) [*c]u8;
pub extern fn a64l(__s: [*c]const u8) c_long;
pub const __u_char = u8;
pub const __u_short = c_ushort;
pub const __u_int = c_uint;
pub const __u_long = c_ulong;
pub const __int8_t = i8;
pub const __uint8_t = u8;
pub const __int16_t = c_short;
pub const __uint16_t = c_ushort;
pub const __int32_t = c_int;
pub const __uint32_t = c_uint;
pub const __int64_t = c_long;
pub const __uint64_t = c_ulong;
pub const __int_least8_t = __int8_t;
pub const __uint_least8_t = __uint8_t;
pub const __int_least16_t = __int16_t;
pub const __uint_least16_t = __uint16_t;
pub const __int_least32_t = __int32_t;
pub const __uint_least32_t = __uint32_t;
pub const __int_least64_t = __int64_t;
pub const __uint_least64_t = __uint64_t;
pub const __quad_t = c_long;
pub const __u_quad_t = c_ulong;
pub const __intmax_t = c_long;
pub const __uintmax_t = c_ulong;
pub const __dev_t = c_ulong;
pub const __uid_t = c_uint;
pub const __gid_t = c_uint;
pub const __ino_t = c_ulong;
pub const __ino64_t = c_ulong;
pub const __mode_t = c_uint;
pub const __nlink_t = c_ulong;
pub const __off_t = c_long;
pub const __off64_t = c_long;
pub const __pid_t = c_int;
pub const __fsid_t = extern struct {
    __val: [2]c_int = @import("std").mem.zeroes([2]c_int),
};
pub const __clock_t = c_long;
pub const __rlim_t = c_ulong;
pub const __rlim64_t = c_ulong;
pub const __id_t = c_uint;
pub const __time_t = c_long;
pub const __useconds_t = c_uint;
pub const __suseconds_t = c_long;
pub const __suseconds64_t = c_long;
pub const __daddr_t = c_int;
pub const __key_t = c_int;
pub const __clockid_t = c_int;
pub const __timer_t = ?*anyopaque;
pub const __blksize_t = c_long;
pub const __blkcnt_t = c_long;
pub const __blkcnt64_t = c_long;
pub const __fsblkcnt_t = c_ulong;
pub const __fsblkcnt64_t = c_ulong;
pub const __fsfilcnt_t = c_ulong;
pub const __fsfilcnt64_t = c_ulong;
pub const __fsword_t = c_long;
pub const __ssize_t = c_long;
pub const __syscall_slong_t = c_long;
pub const __syscall_ulong_t = c_ulong;
pub const __loff_t = __off64_t;
pub const __caddr_t = [*c]u8;
pub const __intptr_t = c_long;
pub const __socklen_t = c_uint;
pub const __sig_atomic_t = c_int;
pub const u_char = __u_char;
pub const u_short = __u_short;
pub const u_int = __u_int;
pub const u_long = __u_long;
pub const quad_t = __quad_t;
pub const u_quad_t = __u_quad_t;
pub const fsid_t = __fsid_t;
pub const loff_t = __loff_t;
pub const ino_t = __ino_t;
pub const dev_t = __dev_t;
pub const gid_t = __gid_t;
pub const mode_t = __mode_t;
pub const nlink_t = __nlink_t;
pub const uid_t = __uid_t;
pub const off_t = __off_t;
pub const pid_t = __pid_t;
pub const id_t = __id_t;
pub const daddr_t = __daddr_t;
pub const caddr_t = __caddr_t;
pub const key_t = __key_t;
pub const clock_t = __clock_t;
pub const clockid_t = __clockid_t;
pub const time_t = __time_t;
pub const timer_t = __timer_t;
pub const ulong = c_ulong;
pub const ushort = c_ushort;
pub const uint = c_uint;
pub const u_int8_t = __uint8_t;
pub const u_int16_t = __uint16_t;
pub const u_int32_t = __uint32_t;
pub const u_int64_t = __uint64_t;
pub const register_t = c_long;
pub fn __bswap_16(arg___bsx: __uint16_t) callconv(.c) __uint16_t {
    var __bsx = arg___bsx;
    _ = &__bsx;
    return @as(__uint16_t, @bitCast(@as(c_short, @truncate(((@as(c_int, @bitCast(@as(c_uint, __bsx))) >> @intCast(8)) & @as(c_int, 255)) | ((@as(c_int, @bitCast(@as(c_uint, __bsx))) & @as(c_int, 255)) << @intCast(8))))));
}
pub fn __bswap_32(arg___bsx: __uint32_t) callconv(.c) __uint32_t {
    var __bsx = arg___bsx;
    _ = &__bsx;
    return ((((__bsx & @as(c_uint, 4278190080)) >> @intCast(24)) | ((__bsx & @as(c_uint, 16711680)) >> @intCast(8))) | ((__bsx & @as(c_uint, 65280)) << @intCast(8))) | ((__bsx & @as(c_uint, 255)) << @intCast(24));
}
pub fn __bswap_64(arg___bsx: __uint64_t) callconv(.c) __uint64_t {
    var __bsx = arg___bsx;
    _ = &__bsx;
    return @as(__uint64_t, @bitCast(@as(c_ulong, @truncate(((((((((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 18374686479671623680)) >> @intCast(56)) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 71776119061217280)) >> @intCast(40))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 280375465082880)) >> @intCast(24))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 1095216660480)) >> @intCast(8))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 4278190080)) << @intCast(8))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 16711680)) << @intCast(24))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 65280)) << @intCast(40))) | ((@as(c_ulonglong, @bitCast(@as(c_ulonglong, __bsx))) & @as(c_ulonglong, 255)) << @intCast(56))))));
}
pub fn __uint16_identity(arg___x: __uint16_t) callconv(.c) __uint16_t {
    var __x = arg___x;
    _ = &__x;
    return __x;
}
pub fn __uint32_identity(arg___x: __uint32_t) callconv(.c) __uint32_t {
    var __x = arg___x;
    _ = &__x;
    return __x;
}
pub fn __uint64_identity(arg___x: __uint64_t) callconv(.c) __uint64_t {
    var __x = arg___x;
    _ = &__x;
    return __x;
}
pub const __sigset_t = extern struct {
    __val: [16]c_ulong = @import("std").mem.zeroes([16]c_ulong),
};
pub const sigset_t = __sigset_t;
pub const struct_timeval = extern struct {
    tv_sec: __time_t = @import("std").mem.zeroes(__time_t),
    tv_usec: __suseconds_t = @import("std").mem.zeroes(__suseconds_t),
};
pub const struct_timespec = extern struct {
    tv_sec: __time_t = @import("std").mem.zeroes(__time_t),
    tv_nsec: __syscall_slong_t = @import("std").mem.zeroes(__syscall_slong_t),
};
pub const suseconds_t = __suseconds_t;
pub const __fd_mask = c_long;
pub const fd_set = extern struct {
    __fds_bits: [16]__fd_mask = @import("std").mem.zeroes([16]__fd_mask),
};
pub const fd_mask = __fd_mask;
pub extern fn select(__nfds: c_int, noalias __readfds: [*c]fd_set, noalias __writefds: [*c]fd_set, noalias __exceptfds: [*c]fd_set, noalias __timeout: [*c]struct_timeval) c_int;
pub extern fn pselect(__nfds: c_int, noalias __readfds: [*c]fd_set, noalias __writefds: [*c]fd_set, noalias __exceptfds: [*c]fd_set, noalias __timeout: [*c]const struct_timespec, noalias __sigmask: [*c]const __sigset_t) c_int;
pub const blksize_t = __blksize_t;
pub const blkcnt_t = __blkcnt_t;
pub const fsblkcnt_t = __fsblkcnt_t;
pub const fsfilcnt_t = __fsfilcnt_t;
const struct_unnamed_1 = extern struct {
    __low: c_uint = @import("std").mem.zeroes(c_uint),
    __high: c_uint = @import("std").mem.zeroes(c_uint),
};
pub const __atomic_wide_counter = extern union {
    __value64: c_ulonglong,
    __value32: struct_unnamed_1,
};
pub const struct___pthread_internal_list = extern struct {
    __prev: [*c]struct___pthread_internal_list = @import("std").mem.zeroes([*c]struct___pthread_internal_list),
    __next: [*c]struct___pthread_internal_list = @import("std").mem.zeroes([*c]struct___pthread_internal_list),
};
pub const __pthread_list_t = struct___pthread_internal_list;
pub const struct___pthread_internal_slist = extern struct {
    __next: [*c]struct___pthread_internal_slist = @import("std").mem.zeroes([*c]struct___pthread_internal_slist),
};
pub const __pthread_slist_t = struct___pthread_internal_slist;
pub const struct___pthread_mutex_s = extern struct {
    __lock: c_int = @import("std").mem.zeroes(c_int),
    __count: c_uint = @import("std").mem.zeroes(c_uint),
    __owner: c_int = @import("std").mem.zeroes(c_int),
    __nusers: c_uint = @import("std").mem.zeroes(c_uint),
    __kind: c_int = @import("std").mem.zeroes(c_int),
    __spins: c_short = @import("std").mem.zeroes(c_short),
    __elision: c_short = @import("std").mem.zeroes(c_short),
    __list: __pthread_list_t = @import("std").mem.zeroes(__pthread_list_t),
};
pub const struct___pthread_rwlock_arch_t = extern struct {
    __readers: c_uint = @import("std").mem.zeroes(c_uint),
    __writers: c_uint = @import("std").mem.zeroes(c_uint),
    __wrphase_futex: c_uint = @import("std").mem.zeroes(c_uint),
    __writers_futex: c_uint = @import("std").mem.zeroes(c_uint),
    __pad3: c_uint = @import("std").mem.zeroes(c_uint),
    __pad4: c_uint = @import("std").mem.zeroes(c_uint),
    __cur_writer: c_int = @import("std").mem.zeroes(c_int),
    __shared: c_int = @import("std").mem.zeroes(c_int),
    __rwelision: i8 = @import("std").mem.zeroes(i8),
    __pad1: [7]u8 = @import("std").mem.zeroes([7]u8),
    __pad2: c_ulong = @import("std").mem.zeroes(c_ulong),
    __flags: c_uint = @import("std").mem.zeroes(c_uint),
};
pub const struct___pthread_cond_s = extern struct {
    __wseq: __atomic_wide_counter = @import("std").mem.zeroes(__atomic_wide_counter),
    __g1_start: __atomic_wide_counter = @import("std").mem.zeroes(__atomic_wide_counter),
    __g_size: [2]c_uint = @import("std").mem.zeroes([2]c_uint),
    __g1_orig_size: c_uint = @import("std").mem.zeroes(c_uint),
    __wrefs: c_uint = @import("std").mem.zeroes(c_uint),
    __g_signals: [2]c_uint = @import("std").mem.zeroes([2]c_uint),
};
pub const __tss_t = c_uint;
pub const __thrd_t = c_ulong;
pub const __once_flag = extern struct {
    __data: c_int = @import("std").mem.zeroes(c_int),
};
pub const pthread_t = c_ulong;
pub const pthread_mutexattr_t = extern union {
    __size: [4]u8,
    __align: c_int,
};
pub const pthread_condattr_t = extern union {
    __size: [4]u8,
    __align: c_int,
};
pub const pthread_key_t = c_uint;
pub const pthread_once_t = c_int;
pub const union_pthread_attr_t = extern union {
    __size: [56]u8,
    __align: c_long,
};
pub const pthread_attr_t = union_pthread_attr_t;
pub const pthread_mutex_t = extern union {
    __data: struct___pthread_mutex_s,
    __size: [40]u8,
    __align: c_long,
};
pub const pthread_cond_t = extern union {
    __data: struct___pthread_cond_s,
    __size: [48]u8,
    __align: c_longlong,
};
pub const pthread_rwlock_t = extern union {
    __data: struct___pthread_rwlock_arch_t,
    __size: [56]u8,
    __align: c_long,
};
pub const pthread_rwlockattr_t = extern union {
    __size: [8]u8,
    __align: c_long,
};
pub const pthread_spinlock_t = c_int;
pub const pthread_barrier_t = extern union {
    __size: [32]u8,
    __align: c_long,
};
pub const pthread_barrierattr_t = extern union {
    __size: [4]u8,
    __align: c_int,
};
pub extern fn random() c_long;
pub extern fn srandom(__seed: c_uint) void;
pub extern fn initstate(__seed: c_uint, __statebuf: [*c]u8, __statelen: usize) [*c]u8;
pub extern fn setstate(__statebuf: [*c]u8) [*c]u8;
pub const struct_random_data = extern struct {
    fptr: [*c]i32 = @import("std").mem.zeroes([*c]i32),
    rptr: [*c]i32 = @import("std").mem.zeroes([*c]i32),
    state: [*c]i32 = @import("std").mem.zeroes([*c]i32),
    rand_type: c_int = @import("std").mem.zeroes(c_int),
    rand_deg: c_int = @import("std").mem.zeroes(c_int),
    rand_sep: c_int = @import("std").mem.zeroes(c_int),
    end_ptr: [*c]i32 = @import("std").mem.zeroes([*c]i32),
};
pub extern fn random_r(noalias __buf: [*c]struct_random_data, noalias __result: [*c]i32) c_int;
pub extern fn srandom_r(__seed: c_uint, __buf: [*c]struct_random_data) c_int;
pub extern fn initstate_r(__seed: c_uint, noalias __statebuf: [*c]u8, __statelen: usize, noalias __buf: [*c]struct_random_data) c_int;
pub extern fn setstate_r(noalias __statebuf: [*c]u8, noalias __buf: [*c]struct_random_data) c_int;
pub extern fn rand() c_int;
pub extern fn srand(__seed: c_uint) void;
pub extern fn rand_r(__seed: [*c]c_uint) c_int;
pub extern fn drand48() f64;
pub extern fn erand48(__xsubi: [*c]c_ushort) f64;
pub extern fn lrand48() c_long;
pub extern fn nrand48(__xsubi: [*c]c_ushort) c_long;
pub extern fn mrand48() c_long;
pub extern fn jrand48(__xsubi: [*c]c_ushort) c_long;
pub extern fn srand48(__seedval: c_long) void;
pub extern fn seed48(__seed16v: [*c]c_ushort) [*c]c_ushort;
pub extern fn lcong48(__param: [*c]c_ushort) void;
pub const struct_drand48_data = extern struct {
    __x: [3]c_ushort = @import("std").mem.zeroes([3]c_ushort),
    __old_x: [3]c_ushort = @import("std").mem.zeroes([3]c_ushort),
    __c: c_ushort = @import("std").mem.zeroes(c_ushort),
    __init: c_ushort = @import("std").mem.zeroes(c_ushort),
    __a: c_ulonglong = @import("std").mem.zeroes(c_ulonglong),
};
pub extern fn drand48_r(noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]f64) c_int;
pub extern fn erand48_r(__xsubi: [*c]c_ushort, noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]f64) c_int;
pub extern fn lrand48_r(noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]c_long) c_int;
pub extern fn nrand48_r(__xsubi: [*c]c_ushort, noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]c_long) c_int;
pub extern fn mrand48_r(noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]c_long) c_int;
pub extern fn jrand48_r(__xsubi: [*c]c_ushort, noalias __buffer: [*c]struct_drand48_data, noalias __result: [*c]c_long) c_int;
pub extern fn srand48_r(__seedval: c_long, __buffer: [*c]struct_drand48_data) c_int;
pub extern fn seed48_r(__seed16v: [*c]c_ushort, __buffer: [*c]struct_drand48_data) c_int;
pub extern fn lcong48_r(__param: [*c]c_ushort, __buffer: [*c]struct_drand48_data) c_int;
pub extern fn arc4random() __uint32_t;
pub extern fn arc4random_buf(__buf: ?*anyopaque, __size: usize) void;
pub extern fn arc4random_uniform(__upper_bound: __uint32_t) __uint32_t;
pub extern fn malloc(__size: c_ulong) ?*anyopaque;
pub extern fn calloc(__nmemb: c_ulong, __size: c_ulong) ?*anyopaque;
pub extern fn realloc(__ptr: ?*anyopaque, __size: c_ulong) ?*anyopaque;
pub extern fn free(__ptr: ?*anyopaque) void;
pub extern fn reallocarray(__ptr: ?*anyopaque, __nmemb: usize, __size: usize) ?*anyopaque;
pub extern fn alloca(__size: c_ulong) ?*anyopaque;
pub extern fn valloc(__size: usize) ?*anyopaque;
pub extern fn posix_memalign(__memptr: [*c]?*anyopaque, __alignment: usize, __size: usize) c_int;
pub extern fn aligned_alloc(__alignment: c_ulong, __size: c_ulong) ?*anyopaque;
pub extern fn abort() noreturn;
pub extern fn atexit(__func: ?*const fn () callconv(.c) void) c_int;
pub extern fn at_quick_exit(__func: ?*const fn () callconv(.c) void) c_int;
pub extern fn on_exit(__func: ?*const fn (c_int, ?*anyopaque) callconv(.c) void, __arg: ?*anyopaque) c_int;
pub extern fn exit(__status: c_int) noreturn;
pub extern fn quick_exit(__status: c_int) noreturn;
pub extern fn _Exit(__status: c_int) noreturn;
pub extern fn getenv(__name: [*c]const u8) [*c]u8;
pub extern fn putenv(__string: [*c]u8) c_int;
pub extern fn setenv(__name: [*c]const u8, __value: [*c]const u8, __replace: c_int) c_int;
pub extern fn unsetenv(__name: [*c]const u8) c_int;
pub extern fn clearenv() c_int;
pub extern fn mktemp(__template: [*c]u8) [*c]u8;
pub extern fn mkstemp(__template: [*c]u8) c_int;
pub extern fn mkstemps(__template: [*c]u8, __suffixlen: c_int) c_int;
pub extern fn mkdtemp(__template: [*c]u8) [*c]u8;
pub extern fn system(__command: [*c]const u8) c_int;
pub extern fn realpath(noalias __name: [*c]const u8, noalias __resolved: [*c]u8) [*c]u8;
pub const __compar_fn_t = ?*const fn (?*const anyopaque, ?*const anyopaque) callconv(.c) c_int;
pub extern fn bsearch(__key: ?*const anyopaque, __base: ?*const anyopaque, __nmemb: usize, __size: usize, __compar: __compar_fn_t) ?*anyopaque;
pub extern fn qsort(__base: ?*anyopaque, __nmemb: usize, __size: usize, __compar: __compar_fn_t) void;
pub extern fn abs(__x: c_int) c_int;
pub extern fn labs(__x: c_long) c_long;
pub extern fn llabs(__x: c_longlong) c_longlong;
pub extern fn div(__numer: c_int, __denom: c_int) div_t;
pub extern fn ldiv(__numer: c_long, __denom: c_long) ldiv_t;
pub extern fn lldiv(__numer: c_longlong, __denom: c_longlong) lldiv_t;
pub extern fn ecvt(__value: f64, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int) [*c]u8;
pub extern fn fcvt(__value: f64, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int) [*c]u8;
pub extern fn gcvt(__value: f64, __ndigit: c_int, __buf: [*c]u8) [*c]u8;
pub extern fn qecvt(__value: c_longdouble, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int) [*c]u8;
pub extern fn qfcvt(__value: c_longdouble, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int) [*c]u8;
pub extern fn qgcvt(__value: c_longdouble, __ndigit: c_int, __buf: [*c]u8) [*c]u8;
pub extern fn ecvt_r(__value: f64, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int, noalias __buf: [*c]u8, __len: usize) c_int;
pub extern fn fcvt_r(__value: f64, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int, noalias __buf: [*c]u8, __len: usize) c_int;
pub extern fn qecvt_r(__value: c_longdouble, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int, noalias __buf: [*c]u8, __len: usize) c_int;
pub extern fn qfcvt_r(__value: c_longdouble, __ndigit: c_int, noalias __decpt: [*c]c_int, noalias __sign: [*c]c_int, noalias __buf: [*c]u8, __len: usize) c_int;
pub extern fn mblen(__s: [*c]const u8, __n: usize) c_int;
pub extern fn mbtowc(noalias __pwc: [*c]wchar_t, noalias __s: [*c]const u8, __n: usize) c_int;
pub extern fn wctomb(__s: [*c]u8, __wchar: wchar_t) c_int;
pub extern fn mbstowcs(noalias __pwcs: [*c]wchar_t, noalias __s: [*c]const u8, __n: usize) usize;
pub extern fn wcstombs(noalias __s: [*c]u8, noalias __pwcs: [*c]const wchar_t, __n: usize) usize;
pub extern fn rpmatch(__response: [*c]const u8) c_int;
pub extern fn getsubopt(noalias __optionp: [*c][*c]u8, noalias __tokens: [*c]const [*c]u8, noalias __valuep: [*c][*c]u8) c_int;
pub extern fn getloadavg(__loadavg: [*c]f64, __nelem: c_int) c_int;
pub const gsl_mode_t = c_uint;
pub const struct___va_list_tag_2 = extern struct {
    gp_offset: c_uint = @import("std").mem.zeroes(c_uint),
    fp_offset: c_uint = @import("std").mem.zeroes(c_uint),
    overflow_arg_area: ?*anyopaque = @import("std").mem.zeroes(?*anyopaque),
    reg_save_area: ?*anyopaque = @import("std").mem.zeroes(?*anyopaque),
};
pub const __builtin_va_list = [1]struct___va_list_tag_2;
pub const __gnuc_va_list = __builtin_va_list;
const union_unnamed_3 = extern union {
    __wch: c_uint,
    __wchb: [4]u8,
};
pub const __mbstate_t = extern struct {
    __count: c_int = @import("std").mem.zeroes(c_int),
    __value: union_unnamed_3 = @import("std").mem.zeroes(union_unnamed_3),
};
pub const struct__G_fpos_t = extern struct {
    __pos: __off_t = @import("std").mem.zeroes(__off_t),
    __state: __mbstate_t = @import("std").mem.zeroes(__mbstate_t),
};
pub const __fpos_t = struct__G_fpos_t;
pub const struct__G_fpos64_t = extern struct {
    __pos: __off64_t = @import("std").mem.zeroes(__off64_t),
    __state: __mbstate_t = @import("std").mem.zeroes(__mbstate_t),
};
pub const __fpos64_t = struct__G_fpos64_t;
pub const struct__IO_marker = opaque {};
// /usr/include/bits/types/struct_FILE.h:74:7: warning: struct demoted to opaque type - has bitfield
pub const struct__IO_FILE = opaque {};
pub const __FILE = struct__IO_FILE;
pub const FILE = struct__IO_FILE;
pub const struct__IO_codecvt = opaque {};
pub const struct__IO_wide_data = opaque {};
pub const _IO_lock_t = anyopaque;
pub const cookie_read_function_t = fn (?*anyopaque, [*c]u8, usize) callconv(.c) __ssize_t;
pub const cookie_write_function_t = fn (?*anyopaque, [*c]const u8, usize) callconv(.c) __ssize_t;
pub const cookie_seek_function_t = fn (?*anyopaque, [*c]__off64_t, c_int) callconv(.c) c_int;
pub const cookie_close_function_t = fn (?*anyopaque) callconv(.c) c_int;
pub const struct__IO_cookie_io_functions_t = extern struct {
    read: ?*const cookie_read_function_t = @import("std").mem.zeroes(?*const cookie_read_function_t),
    write: ?*const cookie_write_function_t = @import("std").mem.zeroes(?*const cookie_write_function_t),
    seek: ?*const cookie_seek_function_t = @import("std").mem.zeroes(?*const cookie_seek_function_t),
    close: ?*const cookie_close_function_t = @import("std").mem.zeroes(?*const cookie_close_function_t),
};
pub const cookie_io_functions_t = struct__IO_cookie_io_functions_t;
pub const va_list = __gnuc_va_list;
pub const fpos_t = __fpos_t;
pub extern var stdin: ?*FILE;
pub extern var stdout: ?*FILE;
pub extern var stderr: ?*FILE;
pub extern fn remove(__filename: [*c]const u8) c_int;
pub extern fn rename(__old: [*c]const u8, __new: [*c]const u8) c_int;
pub extern fn renameat(__oldfd: c_int, __old: [*c]const u8, __newfd: c_int, __new: [*c]const u8) c_int;
pub extern fn fclose(__stream: ?*FILE) c_int;
pub extern fn tmpfile() ?*FILE;
pub extern fn tmpnam([*c]u8) [*c]u8;
pub extern fn tmpnam_r(__s: [*c]u8) [*c]u8;
pub extern fn tempnam(__dir: [*c]const u8, __pfx: [*c]const u8) [*c]u8;
pub extern fn fflush(__stream: ?*FILE) c_int;
pub extern fn fflush_unlocked(__stream: ?*FILE) c_int;
pub extern fn fopen(__filename: [*c]const u8, __modes: [*c]const u8) ?*FILE;
pub extern fn freopen(noalias __filename: [*c]const u8, noalias __modes: [*c]const u8, noalias __stream: ?*FILE) ?*FILE;
pub extern fn fdopen(__fd: c_int, __modes: [*c]const u8) ?*FILE;
pub extern fn fopencookie(noalias __magic_cookie: ?*anyopaque, noalias __modes: [*c]const u8, __io_funcs: cookie_io_functions_t) ?*FILE;
pub extern fn fmemopen(__s: ?*anyopaque, __len: usize, __modes: [*c]const u8) ?*FILE;
pub extern fn open_memstream(__bufloc: [*c][*c]u8, __sizeloc: [*c]usize) ?*FILE;
pub extern fn setbuf(noalias __stream: ?*FILE, noalias __buf: [*c]u8) void;
pub extern fn setvbuf(noalias __stream: ?*FILE, noalias __buf: [*c]u8, __modes: c_int, __n: usize) c_int;
pub extern fn setbuffer(noalias __stream: ?*FILE, noalias __buf: [*c]u8, __size: usize) void;
pub extern fn setlinebuf(__stream: ?*FILE) void;
pub extern fn fprintf(noalias __stream: ?*FILE, noalias __format: [*c]const u8, ...) c_int;
pub extern fn printf(__format: [*c]const u8, ...) c_int;
pub extern fn sprintf(noalias __s: [*c]u8, noalias __format: [*c]const u8, ...) c_int;
pub extern fn vfprintf(noalias __s: ?*FILE, noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn vprintf(noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn vsprintf(noalias __s: [*c]u8, noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn snprintf(noalias __s: [*c]u8, __maxlen: c_ulong, noalias __format: [*c]const u8, ...) c_int;
pub extern fn vsnprintf(noalias __s: [*c]u8, __maxlen: c_ulong, noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn vasprintf(noalias __ptr: [*c][*c]u8, noalias __f: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn __asprintf(noalias __ptr: [*c][*c]u8, noalias __fmt: [*c]const u8, ...) c_int;
pub extern fn asprintf(noalias __ptr: [*c][*c]u8, noalias __fmt: [*c]const u8, ...) c_int;
pub extern fn vdprintf(__fd: c_int, noalias __fmt: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn dprintf(__fd: c_int, noalias __fmt: [*c]const u8, ...) c_int;
pub extern fn fscanf(noalias __stream: ?*FILE, noalias __format: [*c]const u8, ...) c_int;
pub extern fn scanf(noalias __format: [*c]const u8, ...) c_int;
pub extern fn sscanf(noalias __s: [*c]const u8, noalias __format: [*c]const u8, ...) c_int;
pub extern fn vfscanf(noalias __s: ?*FILE, noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn vscanf(noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn vsscanf(noalias __s: [*c]const u8, noalias __format: [*c]const u8, __arg: [*c]struct___va_list_tag_2) c_int;
pub extern fn fgetc(__stream: ?*FILE) c_int;
pub extern fn getc(__stream: ?*FILE) c_int;
pub extern fn getchar() c_int;
pub extern fn getc_unlocked(__stream: ?*FILE) c_int;
pub extern fn getchar_unlocked() c_int;
pub extern fn fgetc_unlocked(__stream: ?*FILE) c_int;
pub extern fn fputc(__c: c_int, __stream: ?*FILE) c_int;
pub extern fn putc(__c: c_int, __stream: ?*FILE) c_int;
pub extern fn putchar(__c: c_int) c_int;
pub extern fn fputc_unlocked(__c: c_int, __stream: ?*FILE) c_int;
pub extern fn putc_unlocked(__c: c_int, __stream: ?*FILE) c_int;
pub extern fn putchar_unlocked(__c: c_int) c_int;
pub extern fn getw(__stream: ?*FILE) c_int;
pub extern fn putw(__w: c_int, __stream: ?*FILE) c_int;
pub extern fn fgets(noalias __s: [*c]u8, __n: c_int, noalias __stream: ?*FILE) [*c]u8;
pub extern fn __getdelim(noalias __lineptr: [*c][*c]u8, noalias __n: [*c]usize, __delimiter: c_int, noalias __stream: ?*FILE) __ssize_t;
pub extern fn getdelim(noalias __lineptr: [*c][*c]u8, noalias __n: [*c]usize, __delimiter: c_int, noalias __stream: ?*FILE) __ssize_t;
pub extern fn getline(noalias __lineptr: [*c][*c]u8, noalias __n: [*c]usize, noalias __stream: ?*FILE) __ssize_t;
pub extern fn fputs(noalias __s: [*c]const u8, noalias __stream: ?*FILE) c_int;
pub extern fn puts(__s: [*c]const u8) c_int;
pub extern fn ungetc(__c: c_int, __stream: ?*FILE) c_int;
pub extern fn fread(__ptr: ?*anyopaque, __size: c_ulong, __n: c_ulong, __stream: ?*FILE) c_ulong;
pub extern fn fwrite(__ptr: ?*const anyopaque, __size: c_ulong, __n: c_ulong, __s: ?*FILE) c_ulong;
pub extern fn fread_unlocked(noalias __ptr: ?*anyopaque, __size: usize, __n: usize, noalias __stream: ?*FILE) usize;
pub extern fn fwrite_unlocked(noalias __ptr: ?*const anyopaque, __size: usize, __n: usize, noalias __stream: ?*FILE) usize;
pub extern fn fseek(__stream: ?*FILE, __off: c_long, __whence: c_int) c_int;
pub extern fn ftell(__stream: ?*FILE) c_long;
pub extern fn rewind(__stream: ?*FILE) void;
pub extern fn fseeko(__stream: ?*FILE, __off: __off_t, __whence: c_int) c_int;
pub extern fn ftello(__stream: ?*FILE) __off_t;
pub extern fn fgetpos(noalias __stream: ?*FILE, noalias __pos: [*c]fpos_t) c_int;
pub extern fn fsetpos(__stream: ?*FILE, __pos: [*c]const fpos_t) c_int;
pub extern fn clearerr(__stream: ?*FILE) void;
pub extern fn feof(__stream: ?*FILE) c_int;
pub extern fn ferror(__stream: ?*FILE) c_int;
pub extern fn clearerr_unlocked(__stream: ?*FILE) void;
pub extern fn feof_unlocked(__stream: ?*FILE) c_int;
pub extern fn ferror_unlocked(__stream: ?*FILE) c_int;
pub extern fn perror(__s: [*c]const u8) void;
pub extern fn fileno(__stream: ?*FILE) c_int;
pub extern fn fileno_unlocked(__stream: ?*FILE) c_int;
pub extern fn pclose(__stream: ?*FILE) c_int;
pub extern fn popen(__command: [*c]const u8, __modes: [*c]const u8) ?*FILE;
pub extern fn ctermid(__s: [*c]u8) [*c]u8;
pub extern fn flockfile(__stream: ?*FILE) void;
pub extern fn ftrylockfile(__stream: ?*FILE) c_int;
pub extern fn funlockfile(__stream: ?*FILE) void;
pub extern fn __uflow(?*FILE) c_int;
pub extern fn __overflow(?*FILE, c_int) c_int;
pub extern fn __errno_location() [*c]c_int;
pub const GSL_SUCCESS: c_int = 0;
pub const GSL_FAILURE: c_int = -1;
pub const GSL_CONTINUE: c_int = -2;
pub const GSL_EDOM: c_int = 1;
pub const GSL_ERANGE: c_int = 2;
pub const GSL_EFAULT: c_int = 3;
pub const GSL_EINVAL: c_int = 4;
pub const GSL_EFAILED: c_int = 5;
pub const GSL_EFACTOR: c_int = 6;
pub const GSL_ESANITY: c_int = 7;
pub const GSL_ENOMEM: c_int = 8;
pub const GSL_EBADFUNC: c_int = 9;
pub const GSL_ERUNAWAY: c_int = 10;
pub const GSL_EMAXITER: c_int = 11;
pub const GSL_EZERODIV: c_int = 12;
pub const GSL_EBADTOL: c_int = 13;
pub const GSL_ETOL: c_int = 14;
pub const GSL_EUNDRFLW: c_int = 15;
pub const GSL_EOVRFLW: c_int = 16;
pub const GSL_ELOSS: c_int = 17;
pub const GSL_EROUND: c_int = 18;
pub const GSL_EBADLEN: c_int = 19;
pub const GSL_ENOTSQR: c_int = 20;
pub const GSL_ESING: c_int = 21;
pub const GSL_EDIVERGE: c_int = 22;
pub const GSL_EUNSUP: c_int = 23;
pub const GSL_EUNIMPL: c_int = 24;
pub const GSL_ECACHE: c_int = 25;
pub const GSL_ETABLE: c_int = 26;
pub const GSL_ENOPROG: c_int = 27;
pub const GSL_ENOPROGJ: c_int = 28;
pub const GSL_ETOLF: c_int = 29;
pub const GSL_ETOLX: c_int = 30;
pub const GSL_ETOLG: c_int = 31;
pub const GSL_EOF: c_int = 32;
const enum_unnamed_4 = c_int;
pub extern fn gsl_error(reason: [*c]const u8, file: [*c]const u8, line: c_int, gsl_errno: c_int) void;
pub extern fn gsl_stream_printf(label: [*c]const u8, file: [*c]const u8, line: c_int, reason: [*c]const u8) void;
pub extern fn gsl_strerror(gsl_errno: c_int) [*c]const u8;
pub const gsl_error_handler_t = fn ([*c]const u8, [*c]const u8, c_int, c_int) callconv(.c) void;
pub const gsl_stream_handler_t = fn ([*c]const u8, [*c]const u8, c_int, [*c]const u8) callconv(.c) void;
pub extern fn gsl_set_error_handler(new_handler: ?*const gsl_error_handler_t) ?*const gsl_error_handler_t;
pub extern fn gsl_set_error_handler_off() ?*const gsl_error_handler_t;
pub extern fn gsl_set_stream_handler(new_handler: ?*const gsl_stream_handler_t) ?*const gsl_stream_handler_t;
pub extern fn gsl_set_stream(new_stream: ?*FILE) ?*FILE;
pub extern var gsl_check_range: c_int;
pub const struct_gsl_permutation_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]usize = @import("std").mem.zeroes([*c]usize),
};
pub const gsl_permutation = struct_gsl_permutation_struct;
pub extern fn gsl_permutation_alloc(n: usize) [*c]gsl_permutation;
pub extern fn gsl_permutation_calloc(n: usize) [*c]gsl_permutation;
pub extern fn gsl_permutation_init(p: [*c]gsl_permutation) void;
pub extern fn gsl_permutation_free(p: [*c]gsl_permutation) void;
pub extern fn gsl_permutation_memcpy(dest: [*c]gsl_permutation, src: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_fread(stream: ?*FILE, p: [*c]gsl_permutation) c_int;
pub extern fn gsl_permutation_fwrite(stream: ?*FILE, p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_fscanf(stream: ?*FILE, p: [*c]gsl_permutation) c_int;
pub extern fn gsl_permutation_fprintf(stream: ?*FILE, p: [*c]const gsl_permutation, format: [*c]const u8) c_int;
pub extern fn gsl_permutation_size(p: [*c]const gsl_permutation) usize;
pub extern fn gsl_permutation_data(p: [*c]const gsl_permutation) [*c]usize;
pub extern fn gsl_permutation_swap(p: [*c]gsl_permutation, i: usize, j: usize) c_int;
pub extern fn gsl_permutation_valid(p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_reverse(p: [*c]gsl_permutation) void;
pub extern fn gsl_permutation_inverse(inv: [*c]gsl_permutation, p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_next(p: [*c]gsl_permutation) c_int;
pub extern fn gsl_permutation_prev(p: [*c]gsl_permutation) c_int;
pub extern fn gsl_permutation_mul(p: [*c]gsl_permutation, pa: [*c]const gsl_permutation, pb: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_linear_to_canonical(q: [*c]gsl_permutation, p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_canonical_to_linear(p: [*c]gsl_permutation, q: [*c]const gsl_permutation) c_int;
pub extern fn gsl_permutation_inversions(p: [*c]const gsl_permutation) usize;
pub extern fn gsl_permutation_linear_cycles(p: [*c]const gsl_permutation) usize;
pub extern fn gsl_permutation_canonical_cycles(q: [*c]const gsl_permutation) usize;
pub extern fn gsl_permutation_get(p: [*c]const gsl_permutation, i: usize) usize;
pub const gsl_complex_packed = [*c]f64;
pub const gsl_complex_packed_float = [*c]f32;
pub const gsl_complex_packed_long_double = [*c]c_longdouble;
pub const gsl_const_complex_packed = [*c]const f64;
pub const gsl_const_complex_packed_float = [*c]const f32;
pub const gsl_const_complex_packed_long_double = [*c]const c_longdouble;
pub const gsl_complex_packed_array = [*c]f64;
pub const gsl_complex_packed_array_float = [*c]f32;
pub const gsl_complex_packed_array_long_double = [*c]c_longdouble;
pub const gsl_const_complex_packed_array = [*c]const f64;
pub const gsl_const_complex_packed_array_float = [*c]const f32;
pub const gsl_const_complex_packed_array_long_double = [*c]const c_longdouble;
pub const gsl_complex_packed_ptr = [*c]f64;
pub const gsl_complex_packed_float_ptr = [*c]f32;
pub const gsl_complex_packed_long_double_ptr = [*c]c_longdouble;
pub const gsl_const_complex_packed_ptr = [*c]const f64;
pub const gsl_const_complex_packed_float_ptr = [*c]const f32;
pub const gsl_const_complex_packed_long_double_ptr = [*c]const c_longdouble;
pub const gsl_complex = extern struct {
    dat: [2]f64 = @import("std").mem.zeroes([2]f64),
};
pub const gsl_complex_long_double = extern struct {
    dat: [2]c_longdouble = @import("std").mem.zeroes([2]c_longdouble),
};
pub const gsl_complex_float = extern struct {
    dat: [2]f32 = @import("std").mem.zeroes([2]f32),
};
pub const struct_gsl_block_long_double_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
};
pub const gsl_block_long_double = struct_gsl_block_long_double_struct;
pub extern fn gsl_block_long_double_alloc(n: usize) [*c]gsl_block_long_double;
pub extern fn gsl_block_long_double_calloc(n: usize) [*c]gsl_block_long_double;
pub extern fn gsl_block_long_double_free(b: [*c]gsl_block_long_double) void;
pub extern fn gsl_block_long_double_fread(stream: ?*FILE, b: [*c]gsl_block_long_double) c_int;
pub extern fn gsl_block_long_double_fwrite(stream: ?*FILE, b: [*c]const gsl_block_long_double) c_int;
pub extern fn gsl_block_long_double_fscanf(stream: ?*FILE, b: [*c]gsl_block_long_double) c_int;
pub extern fn gsl_block_long_double_fprintf(stream: ?*FILE, b: [*c]const gsl_block_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_block_long_double_raw_fread(stream: ?*FILE, b: [*c]c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_double_raw_fwrite(stream: ?*FILE, b: [*c]const c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_double_raw_fscanf(stream: ?*FILE, b: [*c]c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_double_raw_fprintf(stream: ?*FILE, b: [*c]const c_longdouble, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_long_double_size(b: [*c]const gsl_block_long_double) usize;
pub extern fn gsl_block_long_double_data(b: [*c]const gsl_block_long_double) [*c]c_longdouble;
pub const gsl_vector_long_double = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
    block: [*c]gsl_block_long_double = @import("std").mem.zeroes([*c]gsl_block_long_double),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_long_double_view = extern struct {
    vector: gsl_vector_long_double = @import("std").mem.zeroes(gsl_vector_long_double),
};
pub const gsl_vector_long_double_view = _gsl_vector_long_double_view;
pub const _gsl_vector_long_double_const_view = extern struct {
    vector: gsl_vector_long_double = @import("std").mem.zeroes(gsl_vector_long_double),
};
pub const gsl_vector_long_double_const_view = _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_long_double_alloc(n: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_vector_long_double_calloc(n: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_vector_long_double_alloc_from_block(b: [*c]gsl_block_long_double, offset: usize, n: usize, stride: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_vector_long_double_alloc_from_vector(v: [*c]gsl_vector_long_double, offset: usize, n: usize, stride: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_vector_long_double_free(v: [*c]gsl_vector_long_double) void;
pub extern fn gsl_vector_long_double_view_array(v: [*c]c_longdouble, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_vector_long_double_view_array_with_stride(base: [*c]c_longdouble, stride: usize, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_vector_long_double_const_view_array(v: [*c]const c_longdouble, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_long_double_const_view_array_with_stride(base: [*c]const c_longdouble, stride: usize, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_long_double_subvector(v: [*c]gsl_vector_long_double, i: usize, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_vector_long_double_subvector_with_stride(v: [*c]gsl_vector_long_double, i: usize, stride: usize, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_vector_long_double_const_subvector(v: [*c]const gsl_vector_long_double, i: usize, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_long_double_const_subvector_with_stride(v: [*c]const gsl_vector_long_double, i: usize, stride: usize, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_long_double_set_zero(v: [*c]gsl_vector_long_double) void;
pub extern fn gsl_vector_long_double_set_all(v: [*c]gsl_vector_long_double, x: c_longdouble) void;
pub extern fn gsl_vector_long_double_set_basis(v: [*c]gsl_vector_long_double, i: usize) c_int;
pub extern fn gsl_vector_long_double_fread(stream: ?*FILE, v: [*c]gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_fscanf(stream: ?*FILE, v: [*c]gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_vector_long_double_memcpy(dest: [*c]gsl_vector_long_double, src: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_reverse(v: [*c]gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_swap(v: [*c]gsl_vector_long_double, w: [*c]gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_swap_elements(v: [*c]gsl_vector_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_vector_long_double_max(v: [*c]const gsl_vector_long_double) c_longdouble;
pub extern fn gsl_vector_long_double_min(v: [*c]const gsl_vector_long_double) c_longdouble;
pub extern fn gsl_vector_long_double_minmax(v: [*c]const gsl_vector_long_double, min_out: [*c]c_longdouble, max_out: [*c]c_longdouble) void;
pub extern fn gsl_vector_long_double_max_index(v: [*c]const gsl_vector_long_double) usize;
pub extern fn gsl_vector_long_double_min_index(v: [*c]const gsl_vector_long_double) usize;
pub extern fn gsl_vector_long_double_minmax_index(v: [*c]const gsl_vector_long_double, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_long_double_add(a: [*c]gsl_vector_long_double, b: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_sub(a: [*c]gsl_vector_long_double, b: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_mul(a: [*c]gsl_vector_long_double, b: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_div(a: [*c]gsl_vector_long_double, b: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_scale(a: [*c]gsl_vector_long_double, x: c_longdouble) c_int;
pub extern fn gsl_vector_long_double_add_constant(a: [*c]gsl_vector_long_double, x: c_longdouble) c_int;
pub extern fn gsl_vector_long_double_axpby(alpha: c_longdouble, x: [*c]const gsl_vector_long_double, beta: c_longdouble, y: [*c]gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_sum(a: [*c]const gsl_vector_long_double) c_longdouble;
pub extern fn gsl_vector_long_double_equal(u: [*c]const gsl_vector_long_double, v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_isnull(v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_ispos(v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_isneg(v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_isnonneg(v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_long_double_get(v: [*c]const gsl_vector_long_double, i: usize) c_longdouble;
pub extern fn gsl_vector_long_double_set(v: [*c]gsl_vector_long_double, i: usize, x: c_longdouble) void;
pub extern fn gsl_vector_long_double_ptr(v: [*c]gsl_vector_long_double, i: usize) [*c]c_longdouble;
pub extern fn gsl_vector_long_double_const_ptr(v: [*c]const gsl_vector_long_double, i: usize) [*c]const c_longdouble;
pub const struct_gsl_block_complex_long_double_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
};
pub const gsl_block_complex_long_double = struct_gsl_block_complex_long_double_struct;
pub extern fn gsl_block_complex_long_double_alloc(n: usize) [*c]gsl_block_complex_long_double;
pub extern fn gsl_block_complex_long_double_calloc(n: usize) [*c]gsl_block_complex_long_double;
pub extern fn gsl_block_complex_long_double_free(b: [*c]gsl_block_complex_long_double) void;
pub extern fn gsl_block_complex_long_double_fread(stream: ?*FILE, b: [*c]gsl_block_complex_long_double) c_int;
pub extern fn gsl_block_complex_long_double_fwrite(stream: ?*FILE, b: [*c]const gsl_block_complex_long_double) c_int;
pub extern fn gsl_block_complex_long_double_fscanf(stream: ?*FILE, b: [*c]gsl_block_complex_long_double) c_int;
pub extern fn gsl_block_complex_long_double_fprintf(stream: ?*FILE, b: [*c]const gsl_block_complex_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_long_double_raw_fread(stream: ?*FILE, b: [*c]c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_long_double_raw_fwrite(stream: ?*FILE, b: [*c]const c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_long_double_raw_fscanf(stream: ?*FILE, b: [*c]c_longdouble, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_long_double_raw_fprintf(stream: ?*FILE, b: [*c]const c_longdouble, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_long_double_size(b: [*c]const gsl_block_complex_long_double) usize;
pub extern fn gsl_block_complex_long_double_data(b: [*c]const gsl_block_complex_long_double) [*c]c_longdouble;
pub const gsl_vector_complex_long_double = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
    block: [*c]gsl_block_complex_long_double = @import("std").mem.zeroes([*c]gsl_block_complex_long_double),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_complex_long_double_view = extern struct {
    vector: gsl_vector_complex_long_double = @import("std").mem.zeroes(gsl_vector_complex_long_double),
};
pub const gsl_vector_complex_long_double_view = _gsl_vector_complex_long_double_view;
pub const _gsl_vector_complex_long_double_const_view = extern struct {
    vector: gsl_vector_complex_long_double = @import("std").mem.zeroes(gsl_vector_complex_long_double),
};
pub const gsl_vector_complex_long_double_const_view = _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_alloc(n: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_vector_complex_long_double_calloc(n: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_vector_complex_long_double_alloc_from_block(b: [*c]gsl_block_complex_long_double, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_vector_complex_long_double_alloc_from_vector(v: [*c]gsl_vector_complex_long_double, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_vector_complex_long_double_free(v: [*c]gsl_vector_complex_long_double) void;
pub extern fn gsl_vector_complex_long_double_view_array(base: [*c]c_longdouble, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_vector_complex_long_double_view_array_with_stride(base: [*c]c_longdouble, stride: usize, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_vector_complex_long_double_const_view_array(base: [*c]const c_longdouble, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_const_view_array_with_stride(base: [*c]const c_longdouble, stride: usize, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_subvector(base: [*c]gsl_vector_complex_long_double, i: usize, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_vector_complex_long_double_subvector_with_stride(v: [*c]gsl_vector_complex_long_double, i: usize, stride: usize, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_vector_complex_long_double_const_subvector(base: [*c]const gsl_vector_complex_long_double, i: usize, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_const_subvector_with_stride(v: [*c]const gsl_vector_complex_long_double, i: usize, stride: usize, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_real(v: [*c]gsl_vector_complex_long_double) _gsl_vector_long_double_view;
pub extern fn gsl_vector_complex_long_double_imag(v: [*c]gsl_vector_complex_long_double) _gsl_vector_long_double_view;
pub extern fn gsl_vector_complex_long_double_const_real(v: [*c]const gsl_vector_complex_long_double) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_const_imag(v: [*c]const gsl_vector_complex_long_double) _gsl_vector_long_double_const_view;
pub extern fn gsl_vector_complex_long_double_set_zero(v: [*c]gsl_vector_complex_long_double) void;
pub extern fn gsl_vector_complex_long_double_set_all(v: [*c]gsl_vector_complex_long_double, z: gsl_complex_long_double) void;
pub extern fn gsl_vector_complex_long_double_set_basis(v: [*c]gsl_vector_complex_long_double, i: usize) c_int;
pub extern fn gsl_vector_complex_long_double_fread(stream: ?*FILE, v: [*c]gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_fscanf(stream: ?*FILE, v: [*c]gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_complex_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_vector_complex_long_double_memcpy(dest: [*c]gsl_vector_complex_long_double, src: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_conj_memcpy(dest: [*c]gsl_vector_complex_long_double, src: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_reverse(v: [*c]gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_swap(v: [*c]gsl_vector_complex_long_double, w: [*c]gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_swap_elements(v: [*c]gsl_vector_complex_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_vector_complex_long_double_equal(u: [*c]const gsl_vector_complex_long_double, v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_isnull(v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_ispos(v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_isneg(v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_isnonneg(v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_add(a: [*c]gsl_vector_complex_long_double, b: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_sub(a: [*c]gsl_vector_complex_long_double, b: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_mul(a: [*c]gsl_vector_complex_long_double, b: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_div(a: [*c]gsl_vector_complex_long_double, b: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_scale(a: [*c]gsl_vector_complex_long_double, x: gsl_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_add_constant(a: [*c]gsl_vector_complex_long_double, x: gsl_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_axpby(alpha: gsl_complex_long_double, x: [*c]const gsl_vector_complex_long_double, beta: gsl_complex_long_double, y: [*c]gsl_vector_complex_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_div_real(a: [*c]gsl_vector_complex_long_double, b: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_vector_complex_long_double_get(v: [*c]const gsl_vector_complex_long_double, i: usize) gsl_complex_long_double;
pub extern fn gsl_vector_complex_long_double_set(v: [*c]gsl_vector_complex_long_double, i: usize, z: gsl_complex_long_double) void;
pub extern fn gsl_vector_complex_long_double_ptr(v: [*c]gsl_vector_complex_long_double, i: usize) [*c]gsl_complex_long_double;
pub extern fn gsl_vector_complex_long_double_const_ptr(v: [*c]const gsl_vector_complex_long_double, i: usize) [*c]const gsl_complex_long_double;
pub const struct_gsl_block_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
};
pub const gsl_block = struct_gsl_block_struct;
pub extern fn gsl_block_alloc(n: usize) [*c]gsl_block;
pub extern fn gsl_block_calloc(n: usize) [*c]gsl_block;
pub extern fn gsl_block_free(b: [*c]gsl_block) void;
pub extern fn gsl_block_fread(stream: ?*FILE, b: [*c]gsl_block) c_int;
pub extern fn gsl_block_fwrite(stream: ?*FILE, b: [*c]const gsl_block) c_int;
pub extern fn gsl_block_fscanf(stream: ?*FILE, b: [*c]gsl_block) c_int;
pub extern fn gsl_block_fprintf(stream: ?*FILE, b: [*c]const gsl_block, format: [*c]const u8) c_int;
pub extern fn gsl_block_raw_fread(stream: ?*FILE, b: [*c]f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_raw_fwrite(stream: ?*FILE, b: [*c]const f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_raw_fscanf(stream: ?*FILE, b: [*c]f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_raw_fprintf(stream: ?*FILE, b: [*c]const f64, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_size(b: [*c]const gsl_block) usize;
pub extern fn gsl_block_data(b: [*c]const gsl_block) [*c]f64;
pub const gsl_vector = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    block: [*c]gsl_block = @import("std").mem.zeroes([*c]gsl_block),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_view = extern struct {
    vector: gsl_vector = @import("std").mem.zeroes(gsl_vector),
};
pub const gsl_vector_view = _gsl_vector_view;
pub const _gsl_vector_const_view = extern struct {
    vector: gsl_vector = @import("std").mem.zeroes(gsl_vector),
};
pub const gsl_vector_const_view = _gsl_vector_const_view;
pub extern fn gsl_vector_alloc(n: usize) [*c]gsl_vector;
pub extern fn gsl_vector_calloc(n: usize) [*c]gsl_vector;
pub extern fn gsl_vector_alloc_from_block(b: [*c]gsl_block, offset: usize, n: usize, stride: usize) [*c]gsl_vector;
pub extern fn gsl_vector_alloc_from_vector(v: [*c]gsl_vector, offset: usize, n: usize, stride: usize) [*c]gsl_vector;
pub extern fn gsl_vector_free(v: [*c]gsl_vector) void;
pub extern fn gsl_vector_view_array(v: [*c]f64, n: usize) _gsl_vector_view;
pub extern fn gsl_vector_view_array_with_stride(base: [*c]f64, stride: usize, n: usize) _gsl_vector_view;
pub extern fn gsl_vector_const_view_array(v: [*c]const f64, n: usize) _gsl_vector_const_view;
pub extern fn gsl_vector_const_view_array_with_stride(base: [*c]const f64, stride: usize, n: usize) _gsl_vector_const_view;
pub extern fn gsl_vector_subvector(v: [*c]gsl_vector, i: usize, n: usize) _gsl_vector_view;
pub extern fn gsl_vector_subvector_with_stride(v: [*c]gsl_vector, i: usize, stride: usize, n: usize) _gsl_vector_view;
pub extern fn gsl_vector_const_subvector(v: [*c]const gsl_vector, i: usize, n: usize) _gsl_vector_const_view;
pub extern fn gsl_vector_const_subvector_with_stride(v: [*c]const gsl_vector, i: usize, stride: usize, n: usize) _gsl_vector_const_view;
pub extern fn gsl_vector_set_zero(v: [*c]gsl_vector) void;
pub extern fn gsl_vector_set_all(v: [*c]gsl_vector, x: f64) void;
pub extern fn gsl_vector_set_basis(v: [*c]gsl_vector, i: usize) c_int;
pub extern fn gsl_vector_fread(stream: ?*FILE, v: [*c]gsl_vector) c_int;
pub extern fn gsl_vector_fwrite(stream: ?*FILE, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_fscanf(stream: ?*FILE, v: [*c]gsl_vector) c_int;
pub extern fn gsl_vector_fprintf(stream: ?*FILE, v: [*c]const gsl_vector, format: [*c]const u8) c_int;
pub extern fn gsl_vector_memcpy(dest: [*c]gsl_vector, src: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_reverse(v: [*c]gsl_vector) c_int;
pub extern fn gsl_vector_swap(v: [*c]gsl_vector, w: [*c]gsl_vector) c_int;
pub extern fn gsl_vector_swap_elements(v: [*c]gsl_vector, i: usize, j: usize) c_int;
pub extern fn gsl_vector_max(v: [*c]const gsl_vector) f64;
pub extern fn gsl_vector_min(v: [*c]const gsl_vector) f64;
pub extern fn gsl_vector_minmax(v: [*c]const gsl_vector, min_out: [*c]f64, max_out: [*c]f64) void;
pub extern fn gsl_vector_max_index(v: [*c]const gsl_vector) usize;
pub extern fn gsl_vector_min_index(v: [*c]const gsl_vector) usize;
pub extern fn gsl_vector_minmax_index(v: [*c]const gsl_vector, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_add(a: [*c]gsl_vector, b: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_sub(a: [*c]gsl_vector, b: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_mul(a: [*c]gsl_vector, b: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_div(a: [*c]gsl_vector, b: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_scale(a: [*c]gsl_vector, x: f64) c_int;
pub extern fn gsl_vector_add_constant(a: [*c]gsl_vector, x: f64) c_int;
pub extern fn gsl_vector_axpby(alpha: f64, x: [*c]const gsl_vector, beta: f64, y: [*c]gsl_vector) c_int;
pub extern fn gsl_vector_sum(a: [*c]const gsl_vector) f64;
pub extern fn gsl_vector_equal(u: [*c]const gsl_vector, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_isnull(v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_ispos(v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_isneg(v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_isnonneg(v: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_get(v: [*c]const gsl_vector, i: usize) f64;
pub extern fn gsl_vector_set(v: [*c]gsl_vector, i: usize, x: f64) void;
pub extern fn gsl_vector_ptr(v: [*c]gsl_vector, i: usize) [*c]f64;
pub extern fn gsl_vector_const_ptr(v: [*c]const gsl_vector, i: usize) [*c]const f64;
pub const struct_gsl_block_complex_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
};
pub const gsl_block_complex = struct_gsl_block_complex_struct;
pub extern fn gsl_block_complex_alloc(n: usize) [*c]gsl_block_complex;
pub extern fn gsl_block_complex_calloc(n: usize) [*c]gsl_block_complex;
pub extern fn gsl_block_complex_free(b: [*c]gsl_block_complex) void;
pub extern fn gsl_block_complex_fread(stream: ?*FILE, b: [*c]gsl_block_complex) c_int;
pub extern fn gsl_block_complex_fwrite(stream: ?*FILE, b: [*c]const gsl_block_complex) c_int;
pub extern fn gsl_block_complex_fscanf(stream: ?*FILE, b: [*c]gsl_block_complex) c_int;
pub extern fn gsl_block_complex_fprintf(stream: ?*FILE, b: [*c]const gsl_block_complex, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_raw_fread(stream: ?*FILE, b: [*c]f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_raw_fwrite(stream: ?*FILE, b: [*c]const f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_raw_fscanf(stream: ?*FILE, b: [*c]f64, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_raw_fprintf(stream: ?*FILE, b: [*c]const f64, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_size(b: [*c]const gsl_block_complex) usize;
pub extern fn gsl_block_complex_data(b: [*c]const gsl_block_complex) [*c]f64;
pub const gsl_vector_complex = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    block: [*c]gsl_block_complex = @import("std").mem.zeroes([*c]gsl_block_complex),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_complex_view = extern struct {
    vector: gsl_vector_complex = @import("std").mem.zeroes(gsl_vector_complex),
};
pub const gsl_vector_complex_view = _gsl_vector_complex_view;
pub const _gsl_vector_complex_const_view = extern struct {
    vector: gsl_vector_complex = @import("std").mem.zeroes(gsl_vector_complex),
};
pub const gsl_vector_complex_const_view = _gsl_vector_complex_const_view;
pub extern fn gsl_vector_complex_alloc(n: usize) [*c]gsl_vector_complex;
pub extern fn gsl_vector_complex_calloc(n: usize) [*c]gsl_vector_complex;
pub extern fn gsl_vector_complex_alloc_from_block(b: [*c]gsl_block_complex, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex;
pub extern fn gsl_vector_complex_alloc_from_vector(v: [*c]gsl_vector_complex, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex;
pub extern fn gsl_vector_complex_free(v: [*c]gsl_vector_complex) void;
pub extern fn gsl_vector_complex_view_array(base: [*c]f64, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_vector_complex_view_array_with_stride(base: [*c]f64, stride: usize, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_vector_complex_const_view_array(base: [*c]const f64, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_vector_complex_const_view_array_with_stride(base: [*c]const f64, stride: usize, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_vector_complex_subvector(base: [*c]gsl_vector_complex, i: usize, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_vector_complex_subvector_with_stride(v: [*c]gsl_vector_complex, i: usize, stride: usize, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_vector_complex_const_subvector(base: [*c]const gsl_vector_complex, i: usize, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_vector_complex_const_subvector_with_stride(v: [*c]const gsl_vector_complex, i: usize, stride: usize, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_vector_complex_real(v: [*c]gsl_vector_complex) _gsl_vector_view;
pub extern fn gsl_vector_complex_imag(v: [*c]gsl_vector_complex) _gsl_vector_view;
pub extern fn gsl_vector_complex_const_real(v: [*c]const gsl_vector_complex) _gsl_vector_const_view;
pub extern fn gsl_vector_complex_const_imag(v: [*c]const gsl_vector_complex) _gsl_vector_const_view;
pub extern fn gsl_vector_complex_set_zero(v: [*c]gsl_vector_complex) void;
pub extern fn gsl_vector_complex_set_all(v: [*c]gsl_vector_complex, z: gsl_complex) void;
pub extern fn gsl_vector_complex_set_basis(v: [*c]gsl_vector_complex, i: usize) c_int;
pub extern fn gsl_vector_complex_fread(stream: ?*FILE, v: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_fscanf(stream: ?*FILE, v: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_complex, format: [*c]const u8) c_int;
pub extern fn gsl_vector_complex_memcpy(dest: [*c]gsl_vector_complex, src: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_conj_memcpy(dest: [*c]gsl_vector_complex, src: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_reverse(v: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_swap(v: [*c]gsl_vector_complex, w: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_swap_elements(v: [*c]gsl_vector_complex, i: usize, j: usize) c_int;
pub extern fn gsl_vector_complex_equal(u: [*c]const gsl_vector_complex, v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_isnull(v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_ispos(v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_isneg(v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_isnonneg(v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_add(a: [*c]gsl_vector_complex, b: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_sub(a: [*c]gsl_vector_complex, b: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_mul(a: [*c]gsl_vector_complex, b: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_div(a: [*c]gsl_vector_complex, b: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_scale(a: [*c]gsl_vector_complex, x: gsl_complex) c_int;
pub extern fn gsl_vector_complex_add_constant(a: [*c]gsl_vector_complex, x: gsl_complex) c_int;
pub extern fn gsl_vector_complex_axpby(alpha: gsl_complex, x: [*c]const gsl_vector_complex, beta: gsl_complex, y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_vector_complex_div_real(a: [*c]gsl_vector_complex, b: [*c]const gsl_vector) c_int;
pub extern fn gsl_vector_complex_get(v: [*c]const gsl_vector_complex, i: usize) gsl_complex;
pub extern fn gsl_vector_complex_set(v: [*c]gsl_vector_complex, i: usize, z: gsl_complex) void;
pub extern fn gsl_vector_complex_ptr(v: [*c]gsl_vector_complex, i: usize) [*c]gsl_complex;
pub extern fn gsl_vector_complex_const_ptr(v: [*c]const gsl_vector_complex, i: usize) [*c]const gsl_complex;
pub const struct_gsl_block_float_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
};
pub const gsl_block_float = struct_gsl_block_float_struct;
pub extern fn gsl_block_float_alloc(n: usize) [*c]gsl_block_float;
pub extern fn gsl_block_float_calloc(n: usize) [*c]gsl_block_float;
pub extern fn gsl_block_float_free(b: [*c]gsl_block_float) void;
pub extern fn gsl_block_float_fread(stream: ?*FILE, b: [*c]gsl_block_float) c_int;
pub extern fn gsl_block_float_fwrite(stream: ?*FILE, b: [*c]const gsl_block_float) c_int;
pub extern fn gsl_block_float_fscanf(stream: ?*FILE, b: [*c]gsl_block_float) c_int;
pub extern fn gsl_block_float_fprintf(stream: ?*FILE, b: [*c]const gsl_block_float, format: [*c]const u8) c_int;
pub extern fn gsl_block_float_raw_fread(stream: ?*FILE, b: [*c]f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_float_raw_fwrite(stream: ?*FILE, b: [*c]const f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_float_raw_fscanf(stream: ?*FILE, b: [*c]f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_float_raw_fprintf(stream: ?*FILE, b: [*c]const f32, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_float_size(b: [*c]const gsl_block_float) usize;
pub extern fn gsl_block_float_data(b: [*c]const gsl_block_float) [*c]f32;
pub const gsl_vector_float = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
    block: [*c]gsl_block_float = @import("std").mem.zeroes([*c]gsl_block_float),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_float_view = extern struct {
    vector: gsl_vector_float = @import("std").mem.zeroes(gsl_vector_float),
};
pub const gsl_vector_float_view = _gsl_vector_float_view;
pub const _gsl_vector_float_const_view = extern struct {
    vector: gsl_vector_float = @import("std").mem.zeroes(gsl_vector_float),
};
pub const gsl_vector_float_const_view = _gsl_vector_float_const_view;
pub extern fn gsl_vector_float_alloc(n: usize) [*c]gsl_vector_float;
pub extern fn gsl_vector_float_calloc(n: usize) [*c]gsl_vector_float;
pub extern fn gsl_vector_float_alloc_from_block(b: [*c]gsl_block_float, offset: usize, n: usize, stride: usize) [*c]gsl_vector_float;
pub extern fn gsl_vector_float_alloc_from_vector(v: [*c]gsl_vector_float, offset: usize, n: usize, stride: usize) [*c]gsl_vector_float;
pub extern fn gsl_vector_float_free(v: [*c]gsl_vector_float) void;
pub extern fn gsl_vector_float_view_array(v: [*c]f32, n: usize) _gsl_vector_float_view;
pub extern fn gsl_vector_float_view_array_with_stride(base: [*c]f32, stride: usize, n: usize) _gsl_vector_float_view;
pub extern fn gsl_vector_float_const_view_array(v: [*c]const f32, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_vector_float_const_view_array_with_stride(base: [*c]const f32, stride: usize, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_vector_float_subvector(v: [*c]gsl_vector_float, i: usize, n: usize) _gsl_vector_float_view;
pub extern fn gsl_vector_float_subvector_with_stride(v: [*c]gsl_vector_float, i: usize, stride: usize, n: usize) _gsl_vector_float_view;
pub extern fn gsl_vector_float_const_subvector(v: [*c]const gsl_vector_float, i: usize, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_vector_float_const_subvector_with_stride(v: [*c]const gsl_vector_float, i: usize, stride: usize, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_vector_float_set_zero(v: [*c]gsl_vector_float) void;
pub extern fn gsl_vector_float_set_all(v: [*c]gsl_vector_float, x: f32) void;
pub extern fn gsl_vector_float_set_basis(v: [*c]gsl_vector_float, i: usize) c_int;
pub extern fn gsl_vector_float_fread(stream: ?*FILE, v: [*c]gsl_vector_float) c_int;
pub extern fn gsl_vector_float_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_fscanf(stream: ?*FILE, v: [*c]gsl_vector_float) c_int;
pub extern fn gsl_vector_float_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_float, format: [*c]const u8) c_int;
pub extern fn gsl_vector_float_memcpy(dest: [*c]gsl_vector_float, src: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_reverse(v: [*c]gsl_vector_float) c_int;
pub extern fn gsl_vector_float_swap(v: [*c]gsl_vector_float, w: [*c]gsl_vector_float) c_int;
pub extern fn gsl_vector_float_swap_elements(v: [*c]gsl_vector_float, i: usize, j: usize) c_int;
pub extern fn gsl_vector_float_max(v: [*c]const gsl_vector_float) f32;
pub extern fn gsl_vector_float_min(v: [*c]const gsl_vector_float) f32;
pub extern fn gsl_vector_float_minmax(v: [*c]const gsl_vector_float, min_out: [*c]f32, max_out: [*c]f32) void;
pub extern fn gsl_vector_float_max_index(v: [*c]const gsl_vector_float) usize;
pub extern fn gsl_vector_float_min_index(v: [*c]const gsl_vector_float) usize;
pub extern fn gsl_vector_float_minmax_index(v: [*c]const gsl_vector_float, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_float_add(a: [*c]gsl_vector_float, b: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_sub(a: [*c]gsl_vector_float, b: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_mul(a: [*c]gsl_vector_float, b: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_div(a: [*c]gsl_vector_float, b: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_scale(a: [*c]gsl_vector_float, x: f32) c_int;
pub extern fn gsl_vector_float_add_constant(a: [*c]gsl_vector_float, x: f32) c_int;
pub extern fn gsl_vector_float_axpby(alpha: f32, x: [*c]const gsl_vector_float, beta: f32, y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_vector_float_sum(a: [*c]const gsl_vector_float) f32;
pub extern fn gsl_vector_float_equal(u: [*c]const gsl_vector_float, v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_isnull(v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_ispos(v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_isneg(v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_isnonneg(v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_float_get(v: [*c]const gsl_vector_float, i: usize) f32;
pub extern fn gsl_vector_float_set(v: [*c]gsl_vector_float, i: usize, x: f32) void;
pub extern fn gsl_vector_float_ptr(v: [*c]gsl_vector_float, i: usize) [*c]f32;
pub extern fn gsl_vector_float_const_ptr(v: [*c]const gsl_vector_float, i: usize) [*c]const f32;
pub const struct_gsl_block_complex_float_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
};
pub const gsl_block_complex_float = struct_gsl_block_complex_float_struct;
pub extern fn gsl_block_complex_float_alloc(n: usize) [*c]gsl_block_complex_float;
pub extern fn gsl_block_complex_float_calloc(n: usize) [*c]gsl_block_complex_float;
pub extern fn gsl_block_complex_float_free(b: [*c]gsl_block_complex_float) void;
pub extern fn gsl_block_complex_float_fread(stream: ?*FILE, b: [*c]gsl_block_complex_float) c_int;
pub extern fn gsl_block_complex_float_fwrite(stream: ?*FILE, b: [*c]const gsl_block_complex_float) c_int;
pub extern fn gsl_block_complex_float_fscanf(stream: ?*FILE, b: [*c]gsl_block_complex_float) c_int;
pub extern fn gsl_block_complex_float_fprintf(stream: ?*FILE, b: [*c]const gsl_block_complex_float, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_float_raw_fread(stream: ?*FILE, b: [*c]f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_float_raw_fwrite(stream: ?*FILE, b: [*c]const f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_float_raw_fscanf(stream: ?*FILE, b: [*c]f32, n: usize, stride: usize) c_int;
pub extern fn gsl_block_complex_float_raw_fprintf(stream: ?*FILE, b: [*c]const f32, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_complex_float_size(b: [*c]const gsl_block_complex_float) usize;
pub extern fn gsl_block_complex_float_data(b: [*c]const gsl_block_complex_float) [*c]f32;
pub const gsl_vector_complex_float = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
    block: [*c]gsl_block_complex_float = @import("std").mem.zeroes([*c]gsl_block_complex_float),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_complex_float_view = extern struct {
    vector: gsl_vector_complex_float = @import("std").mem.zeroes(gsl_vector_complex_float),
};
pub const gsl_vector_complex_float_view = _gsl_vector_complex_float_view;
pub const _gsl_vector_complex_float_const_view = extern struct {
    vector: gsl_vector_complex_float = @import("std").mem.zeroes(gsl_vector_complex_float),
};
pub const gsl_vector_complex_float_const_view = _gsl_vector_complex_float_const_view;
pub extern fn gsl_vector_complex_float_alloc(n: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_vector_complex_float_calloc(n: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_vector_complex_float_alloc_from_block(b: [*c]gsl_block_complex_float, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_vector_complex_float_alloc_from_vector(v: [*c]gsl_vector_complex_float, offset: usize, n: usize, stride: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_vector_complex_float_free(v: [*c]gsl_vector_complex_float) void;
pub extern fn gsl_vector_complex_float_view_array(base: [*c]f32, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_vector_complex_float_view_array_with_stride(base: [*c]f32, stride: usize, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_vector_complex_float_const_view_array(base: [*c]const f32, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_vector_complex_float_const_view_array_with_stride(base: [*c]const f32, stride: usize, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_vector_complex_float_subvector(base: [*c]gsl_vector_complex_float, i: usize, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_vector_complex_float_subvector_with_stride(v: [*c]gsl_vector_complex_float, i: usize, stride: usize, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_vector_complex_float_const_subvector(base: [*c]const gsl_vector_complex_float, i: usize, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_vector_complex_float_const_subvector_with_stride(v: [*c]const gsl_vector_complex_float, i: usize, stride: usize, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_vector_complex_float_real(v: [*c]gsl_vector_complex_float) _gsl_vector_float_view;
pub extern fn gsl_vector_complex_float_imag(v: [*c]gsl_vector_complex_float) _gsl_vector_float_view;
pub extern fn gsl_vector_complex_float_const_real(v: [*c]const gsl_vector_complex_float) _gsl_vector_float_const_view;
pub extern fn gsl_vector_complex_float_const_imag(v: [*c]const gsl_vector_complex_float) _gsl_vector_float_const_view;
pub extern fn gsl_vector_complex_float_set_zero(v: [*c]gsl_vector_complex_float) void;
pub extern fn gsl_vector_complex_float_set_all(v: [*c]gsl_vector_complex_float, z: gsl_complex_float) void;
pub extern fn gsl_vector_complex_float_set_basis(v: [*c]gsl_vector_complex_float, i: usize) c_int;
pub extern fn gsl_vector_complex_float_fread(stream: ?*FILE, v: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_fscanf(stream: ?*FILE, v: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_complex_float, format: [*c]const u8) c_int;
pub extern fn gsl_vector_complex_float_memcpy(dest: [*c]gsl_vector_complex_float, src: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_conj_memcpy(dest: [*c]gsl_vector_complex_float, src: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_reverse(v: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_swap(v: [*c]gsl_vector_complex_float, w: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_swap_elements(v: [*c]gsl_vector_complex_float, i: usize, j: usize) c_int;
pub extern fn gsl_vector_complex_float_equal(u: [*c]const gsl_vector_complex_float, v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_isnull(v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_ispos(v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_isneg(v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_isnonneg(v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_add(a: [*c]gsl_vector_complex_float, b: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_sub(a: [*c]gsl_vector_complex_float, b: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_mul(a: [*c]gsl_vector_complex_float, b: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_div(a: [*c]gsl_vector_complex_float, b: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_scale(a: [*c]gsl_vector_complex_float, x: gsl_complex_float) c_int;
pub extern fn gsl_vector_complex_float_add_constant(a: [*c]gsl_vector_complex_float, x: gsl_complex_float) c_int;
pub extern fn gsl_vector_complex_float_axpby(alpha: gsl_complex_float, x: [*c]const gsl_vector_complex_float, beta: gsl_complex_float, y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_vector_complex_float_div_real(a: [*c]gsl_vector_complex_float, b: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_vector_complex_float_get(v: [*c]const gsl_vector_complex_float, i: usize) gsl_complex_float;
pub extern fn gsl_vector_complex_float_set(v: [*c]gsl_vector_complex_float, i: usize, z: gsl_complex_float) void;
pub extern fn gsl_vector_complex_float_ptr(v: [*c]gsl_vector_complex_float, i: usize) [*c]gsl_complex_float;
pub extern fn gsl_vector_complex_float_const_ptr(v: [*c]const gsl_vector_complex_float, i: usize) [*c]const gsl_complex_float;
pub const struct_gsl_block_ulong_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ulong = @import("std").mem.zeroes([*c]c_ulong),
};
pub const gsl_block_ulong = struct_gsl_block_ulong_struct;
pub extern fn gsl_block_ulong_alloc(n: usize) [*c]gsl_block_ulong;
pub extern fn gsl_block_ulong_calloc(n: usize) [*c]gsl_block_ulong;
pub extern fn gsl_block_ulong_free(b: [*c]gsl_block_ulong) void;
pub extern fn gsl_block_ulong_fread(stream: ?*FILE, b: [*c]gsl_block_ulong) c_int;
pub extern fn gsl_block_ulong_fwrite(stream: ?*FILE, b: [*c]const gsl_block_ulong) c_int;
pub extern fn gsl_block_ulong_fscanf(stream: ?*FILE, b: [*c]gsl_block_ulong) c_int;
pub extern fn gsl_block_ulong_fprintf(stream: ?*FILE, b: [*c]const gsl_block_ulong, format: [*c]const u8) c_int;
pub extern fn gsl_block_ulong_raw_fread(stream: ?*FILE, b: [*c]c_ulong, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ulong_raw_fwrite(stream: ?*FILE, b: [*c]const c_ulong, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ulong_raw_fscanf(stream: ?*FILE, b: [*c]c_ulong, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ulong_raw_fprintf(stream: ?*FILE, b: [*c]const c_ulong, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_ulong_size(b: [*c]const gsl_block_ulong) usize;
pub extern fn gsl_block_ulong_data(b: [*c]const gsl_block_ulong) [*c]c_ulong;
pub const gsl_vector_ulong = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ulong = @import("std").mem.zeroes([*c]c_ulong),
    block: [*c]gsl_block_ulong = @import("std").mem.zeroes([*c]gsl_block_ulong),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_ulong_view = extern struct {
    vector: gsl_vector_ulong = @import("std").mem.zeroes(gsl_vector_ulong),
};
pub const gsl_vector_ulong_view = _gsl_vector_ulong_view;
pub const _gsl_vector_ulong_const_view = extern struct {
    vector: gsl_vector_ulong = @import("std").mem.zeroes(gsl_vector_ulong),
};
pub const gsl_vector_ulong_const_view = _gsl_vector_ulong_const_view;
pub extern fn gsl_vector_ulong_alloc(n: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_vector_ulong_calloc(n: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_vector_ulong_alloc_from_block(b: [*c]gsl_block_ulong, offset: usize, n: usize, stride: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_vector_ulong_alloc_from_vector(v: [*c]gsl_vector_ulong, offset: usize, n: usize, stride: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_vector_ulong_free(v: [*c]gsl_vector_ulong) void;
pub extern fn gsl_vector_ulong_view_array(v: [*c]c_ulong, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_vector_ulong_view_array_with_stride(base: [*c]c_ulong, stride: usize, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_vector_ulong_const_view_array(v: [*c]const c_ulong, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_vector_ulong_const_view_array_with_stride(base: [*c]const c_ulong, stride: usize, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_vector_ulong_subvector(v: [*c]gsl_vector_ulong, i: usize, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_vector_ulong_subvector_with_stride(v: [*c]gsl_vector_ulong, i: usize, stride: usize, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_vector_ulong_const_subvector(v: [*c]const gsl_vector_ulong, i: usize, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_vector_ulong_const_subvector_with_stride(v: [*c]const gsl_vector_ulong, i: usize, stride: usize, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_vector_ulong_set_zero(v: [*c]gsl_vector_ulong) void;
pub extern fn gsl_vector_ulong_set_all(v: [*c]gsl_vector_ulong, x: c_ulong) void;
pub extern fn gsl_vector_ulong_set_basis(v: [*c]gsl_vector_ulong, i: usize) c_int;
pub extern fn gsl_vector_ulong_fread(stream: ?*FILE, v: [*c]gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_fscanf(stream: ?*FILE, v: [*c]gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_ulong, format: [*c]const u8) c_int;
pub extern fn gsl_vector_ulong_memcpy(dest: [*c]gsl_vector_ulong, src: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_reverse(v: [*c]gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_swap(v: [*c]gsl_vector_ulong, w: [*c]gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_swap_elements(v: [*c]gsl_vector_ulong, i: usize, j: usize) c_int;
pub extern fn gsl_vector_ulong_max(v: [*c]const gsl_vector_ulong) c_ulong;
pub extern fn gsl_vector_ulong_min(v: [*c]const gsl_vector_ulong) c_ulong;
pub extern fn gsl_vector_ulong_minmax(v: [*c]const gsl_vector_ulong, min_out: [*c]c_ulong, max_out: [*c]c_ulong) void;
pub extern fn gsl_vector_ulong_max_index(v: [*c]const gsl_vector_ulong) usize;
pub extern fn gsl_vector_ulong_min_index(v: [*c]const gsl_vector_ulong) usize;
pub extern fn gsl_vector_ulong_minmax_index(v: [*c]const gsl_vector_ulong, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_ulong_add(a: [*c]gsl_vector_ulong, b: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_sub(a: [*c]gsl_vector_ulong, b: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_mul(a: [*c]gsl_vector_ulong, b: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_div(a: [*c]gsl_vector_ulong, b: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_scale(a: [*c]gsl_vector_ulong, x: c_ulong) c_int;
pub extern fn gsl_vector_ulong_add_constant(a: [*c]gsl_vector_ulong, x: c_ulong) c_int;
pub extern fn gsl_vector_ulong_axpby(alpha: c_ulong, x: [*c]const gsl_vector_ulong, beta: c_ulong, y: [*c]gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_sum(a: [*c]const gsl_vector_ulong) c_ulong;
pub extern fn gsl_vector_ulong_equal(u: [*c]const gsl_vector_ulong, v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_isnull(v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_ispos(v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_isneg(v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_isnonneg(v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_vector_ulong_get(v: [*c]const gsl_vector_ulong, i: usize) c_ulong;
pub extern fn gsl_vector_ulong_set(v: [*c]gsl_vector_ulong, i: usize, x: c_ulong) void;
pub extern fn gsl_vector_ulong_ptr(v: [*c]gsl_vector_ulong, i: usize) [*c]c_ulong;
pub extern fn gsl_vector_ulong_const_ptr(v: [*c]const gsl_vector_ulong, i: usize) [*c]const c_ulong;
pub const struct_gsl_block_long_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_long = @import("std").mem.zeroes([*c]c_long),
};
pub const gsl_block_long = struct_gsl_block_long_struct;
pub extern fn gsl_block_long_alloc(n: usize) [*c]gsl_block_long;
pub extern fn gsl_block_long_calloc(n: usize) [*c]gsl_block_long;
pub extern fn gsl_block_long_free(b: [*c]gsl_block_long) void;
pub extern fn gsl_block_long_fread(stream: ?*FILE, b: [*c]gsl_block_long) c_int;
pub extern fn gsl_block_long_fwrite(stream: ?*FILE, b: [*c]const gsl_block_long) c_int;
pub extern fn gsl_block_long_fscanf(stream: ?*FILE, b: [*c]gsl_block_long) c_int;
pub extern fn gsl_block_long_fprintf(stream: ?*FILE, b: [*c]const gsl_block_long, format: [*c]const u8) c_int;
pub extern fn gsl_block_long_raw_fread(stream: ?*FILE, b: [*c]c_long, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_raw_fwrite(stream: ?*FILE, b: [*c]const c_long, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_raw_fscanf(stream: ?*FILE, b: [*c]c_long, n: usize, stride: usize) c_int;
pub extern fn gsl_block_long_raw_fprintf(stream: ?*FILE, b: [*c]const c_long, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_long_size(b: [*c]const gsl_block_long) usize;
pub extern fn gsl_block_long_data(b: [*c]const gsl_block_long) [*c]c_long;
pub const gsl_vector_long = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_long = @import("std").mem.zeroes([*c]c_long),
    block: [*c]gsl_block_long = @import("std").mem.zeroes([*c]gsl_block_long),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_long_view = extern struct {
    vector: gsl_vector_long = @import("std").mem.zeroes(gsl_vector_long),
};
pub const gsl_vector_long_view = _gsl_vector_long_view;
pub const _gsl_vector_long_const_view = extern struct {
    vector: gsl_vector_long = @import("std").mem.zeroes(gsl_vector_long),
};
pub const gsl_vector_long_const_view = _gsl_vector_long_const_view;
pub extern fn gsl_vector_long_alloc(n: usize) [*c]gsl_vector_long;
pub extern fn gsl_vector_long_calloc(n: usize) [*c]gsl_vector_long;
pub extern fn gsl_vector_long_alloc_from_block(b: [*c]gsl_block_long, offset: usize, n: usize, stride: usize) [*c]gsl_vector_long;
pub extern fn gsl_vector_long_alloc_from_vector(v: [*c]gsl_vector_long, offset: usize, n: usize, stride: usize) [*c]gsl_vector_long;
pub extern fn gsl_vector_long_free(v: [*c]gsl_vector_long) void;
pub extern fn gsl_vector_long_view_array(v: [*c]c_long, n: usize) _gsl_vector_long_view;
pub extern fn gsl_vector_long_view_array_with_stride(base: [*c]c_long, stride: usize, n: usize) _gsl_vector_long_view;
pub extern fn gsl_vector_long_const_view_array(v: [*c]const c_long, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_vector_long_const_view_array_with_stride(base: [*c]const c_long, stride: usize, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_vector_long_subvector(v: [*c]gsl_vector_long, i: usize, n: usize) _gsl_vector_long_view;
pub extern fn gsl_vector_long_subvector_with_stride(v: [*c]gsl_vector_long, i: usize, stride: usize, n: usize) _gsl_vector_long_view;
pub extern fn gsl_vector_long_const_subvector(v: [*c]const gsl_vector_long, i: usize, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_vector_long_const_subvector_with_stride(v: [*c]const gsl_vector_long, i: usize, stride: usize, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_vector_long_set_zero(v: [*c]gsl_vector_long) void;
pub extern fn gsl_vector_long_set_all(v: [*c]gsl_vector_long, x: c_long) void;
pub extern fn gsl_vector_long_set_basis(v: [*c]gsl_vector_long, i: usize) c_int;
pub extern fn gsl_vector_long_fread(stream: ?*FILE, v: [*c]gsl_vector_long) c_int;
pub extern fn gsl_vector_long_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_fscanf(stream: ?*FILE, v: [*c]gsl_vector_long) c_int;
pub extern fn gsl_vector_long_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_long, format: [*c]const u8) c_int;
pub extern fn gsl_vector_long_memcpy(dest: [*c]gsl_vector_long, src: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_reverse(v: [*c]gsl_vector_long) c_int;
pub extern fn gsl_vector_long_swap(v: [*c]gsl_vector_long, w: [*c]gsl_vector_long) c_int;
pub extern fn gsl_vector_long_swap_elements(v: [*c]gsl_vector_long, i: usize, j: usize) c_int;
pub extern fn gsl_vector_long_max(v: [*c]const gsl_vector_long) c_long;
pub extern fn gsl_vector_long_min(v: [*c]const gsl_vector_long) c_long;
pub extern fn gsl_vector_long_minmax(v: [*c]const gsl_vector_long, min_out: [*c]c_long, max_out: [*c]c_long) void;
pub extern fn gsl_vector_long_max_index(v: [*c]const gsl_vector_long) usize;
pub extern fn gsl_vector_long_min_index(v: [*c]const gsl_vector_long) usize;
pub extern fn gsl_vector_long_minmax_index(v: [*c]const gsl_vector_long, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_long_add(a: [*c]gsl_vector_long, b: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_sub(a: [*c]gsl_vector_long, b: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_mul(a: [*c]gsl_vector_long, b: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_div(a: [*c]gsl_vector_long, b: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_scale(a: [*c]gsl_vector_long, x: c_long) c_int;
pub extern fn gsl_vector_long_add_constant(a: [*c]gsl_vector_long, x: c_long) c_int;
pub extern fn gsl_vector_long_axpby(alpha: c_long, x: [*c]const gsl_vector_long, beta: c_long, y: [*c]gsl_vector_long) c_int;
pub extern fn gsl_vector_long_sum(a: [*c]const gsl_vector_long) c_long;
pub extern fn gsl_vector_long_equal(u: [*c]const gsl_vector_long, v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_isnull(v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_ispos(v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_isneg(v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_isnonneg(v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_vector_long_get(v: [*c]const gsl_vector_long, i: usize) c_long;
pub extern fn gsl_vector_long_set(v: [*c]gsl_vector_long, i: usize, x: c_long) void;
pub extern fn gsl_vector_long_ptr(v: [*c]gsl_vector_long, i: usize) [*c]c_long;
pub extern fn gsl_vector_long_const_ptr(v: [*c]const gsl_vector_long, i: usize) [*c]const c_long;
pub const struct_gsl_block_uint_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_uint = @import("std").mem.zeroes([*c]c_uint),
};
pub const gsl_block_uint = struct_gsl_block_uint_struct;
pub extern fn gsl_block_uint_alloc(n: usize) [*c]gsl_block_uint;
pub extern fn gsl_block_uint_calloc(n: usize) [*c]gsl_block_uint;
pub extern fn gsl_block_uint_free(b: [*c]gsl_block_uint) void;
pub extern fn gsl_block_uint_fread(stream: ?*FILE, b: [*c]gsl_block_uint) c_int;
pub extern fn gsl_block_uint_fwrite(stream: ?*FILE, b: [*c]const gsl_block_uint) c_int;
pub extern fn gsl_block_uint_fscanf(stream: ?*FILE, b: [*c]gsl_block_uint) c_int;
pub extern fn gsl_block_uint_fprintf(stream: ?*FILE, b: [*c]const gsl_block_uint, format: [*c]const u8) c_int;
pub extern fn gsl_block_uint_raw_fread(stream: ?*FILE, b: [*c]c_uint, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uint_raw_fwrite(stream: ?*FILE, b: [*c]const c_uint, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uint_raw_fscanf(stream: ?*FILE, b: [*c]c_uint, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uint_raw_fprintf(stream: ?*FILE, b: [*c]const c_uint, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_uint_size(b: [*c]const gsl_block_uint) usize;
pub extern fn gsl_block_uint_data(b: [*c]const gsl_block_uint) [*c]c_uint;
pub const gsl_vector_uint = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_uint = @import("std").mem.zeroes([*c]c_uint),
    block: [*c]gsl_block_uint = @import("std").mem.zeroes([*c]gsl_block_uint),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_uint_view = extern struct {
    vector: gsl_vector_uint = @import("std").mem.zeroes(gsl_vector_uint),
};
pub const gsl_vector_uint_view = _gsl_vector_uint_view;
pub const _gsl_vector_uint_const_view = extern struct {
    vector: gsl_vector_uint = @import("std").mem.zeroes(gsl_vector_uint),
};
pub const gsl_vector_uint_const_view = _gsl_vector_uint_const_view;
pub extern fn gsl_vector_uint_alloc(n: usize) [*c]gsl_vector_uint;
pub extern fn gsl_vector_uint_calloc(n: usize) [*c]gsl_vector_uint;
pub extern fn gsl_vector_uint_alloc_from_block(b: [*c]gsl_block_uint, offset: usize, n: usize, stride: usize) [*c]gsl_vector_uint;
pub extern fn gsl_vector_uint_alloc_from_vector(v: [*c]gsl_vector_uint, offset: usize, n: usize, stride: usize) [*c]gsl_vector_uint;
pub extern fn gsl_vector_uint_free(v: [*c]gsl_vector_uint) void;
pub extern fn gsl_vector_uint_view_array(v: [*c]c_uint, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_vector_uint_view_array_with_stride(base: [*c]c_uint, stride: usize, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_vector_uint_const_view_array(v: [*c]const c_uint, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_vector_uint_const_view_array_with_stride(base: [*c]const c_uint, stride: usize, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_vector_uint_subvector(v: [*c]gsl_vector_uint, i: usize, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_vector_uint_subvector_with_stride(v: [*c]gsl_vector_uint, i: usize, stride: usize, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_vector_uint_const_subvector(v: [*c]const gsl_vector_uint, i: usize, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_vector_uint_const_subvector_with_stride(v: [*c]const gsl_vector_uint, i: usize, stride: usize, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_vector_uint_set_zero(v: [*c]gsl_vector_uint) void;
pub extern fn gsl_vector_uint_set_all(v: [*c]gsl_vector_uint, x: c_uint) void;
pub extern fn gsl_vector_uint_set_basis(v: [*c]gsl_vector_uint, i: usize) c_int;
pub extern fn gsl_vector_uint_fread(stream: ?*FILE, v: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_fscanf(stream: ?*FILE, v: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_uint, format: [*c]const u8) c_int;
pub extern fn gsl_vector_uint_memcpy(dest: [*c]gsl_vector_uint, src: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_reverse(v: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_swap(v: [*c]gsl_vector_uint, w: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_swap_elements(v: [*c]gsl_vector_uint, i: usize, j: usize) c_int;
pub extern fn gsl_vector_uint_max(v: [*c]const gsl_vector_uint) c_uint;
pub extern fn gsl_vector_uint_min(v: [*c]const gsl_vector_uint) c_uint;
pub extern fn gsl_vector_uint_minmax(v: [*c]const gsl_vector_uint, min_out: [*c]c_uint, max_out: [*c]c_uint) void;
pub extern fn gsl_vector_uint_max_index(v: [*c]const gsl_vector_uint) usize;
pub extern fn gsl_vector_uint_min_index(v: [*c]const gsl_vector_uint) usize;
pub extern fn gsl_vector_uint_minmax_index(v: [*c]const gsl_vector_uint, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_uint_add(a: [*c]gsl_vector_uint, b: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_sub(a: [*c]gsl_vector_uint, b: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_mul(a: [*c]gsl_vector_uint, b: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_div(a: [*c]gsl_vector_uint, b: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_scale(a: [*c]gsl_vector_uint, x: c_uint) c_int;
pub extern fn gsl_vector_uint_add_constant(a: [*c]gsl_vector_uint, x: c_uint) c_int;
pub extern fn gsl_vector_uint_axpby(alpha: c_uint, x: [*c]const gsl_vector_uint, beta: c_uint, y: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_sum(a: [*c]const gsl_vector_uint) c_uint;
pub extern fn gsl_vector_uint_equal(u: [*c]const gsl_vector_uint, v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_isnull(v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_ispos(v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_isneg(v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_isnonneg(v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_vector_uint_get(v: [*c]const gsl_vector_uint, i: usize) c_uint;
pub extern fn gsl_vector_uint_set(v: [*c]gsl_vector_uint, i: usize, x: c_uint) void;
pub extern fn gsl_vector_uint_ptr(v: [*c]gsl_vector_uint, i: usize) [*c]c_uint;
pub extern fn gsl_vector_uint_const_ptr(v: [*c]const gsl_vector_uint, i: usize) [*c]const c_uint;
pub const struct_gsl_block_int_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_int = @import("std").mem.zeroes([*c]c_int),
};
pub const gsl_block_int = struct_gsl_block_int_struct;
pub extern fn gsl_block_int_alloc(n: usize) [*c]gsl_block_int;
pub extern fn gsl_block_int_calloc(n: usize) [*c]gsl_block_int;
pub extern fn gsl_block_int_free(b: [*c]gsl_block_int) void;
pub extern fn gsl_block_int_fread(stream: ?*FILE, b: [*c]gsl_block_int) c_int;
pub extern fn gsl_block_int_fwrite(stream: ?*FILE, b: [*c]const gsl_block_int) c_int;
pub extern fn gsl_block_int_fscanf(stream: ?*FILE, b: [*c]gsl_block_int) c_int;
pub extern fn gsl_block_int_fprintf(stream: ?*FILE, b: [*c]const gsl_block_int, format: [*c]const u8) c_int;
pub extern fn gsl_block_int_raw_fread(stream: ?*FILE, b: [*c]c_int, n: usize, stride: usize) c_int;
pub extern fn gsl_block_int_raw_fwrite(stream: ?*FILE, b: [*c]const c_int, n: usize, stride: usize) c_int;
pub extern fn gsl_block_int_raw_fscanf(stream: ?*FILE, b: [*c]c_int, n: usize, stride: usize) c_int;
pub extern fn gsl_block_int_raw_fprintf(stream: ?*FILE, b: [*c]const c_int, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_int_size(b: [*c]const gsl_block_int) usize;
pub extern fn gsl_block_int_data(b: [*c]const gsl_block_int) [*c]c_int;
pub const gsl_vector_int = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_int = @import("std").mem.zeroes([*c]c_int),
    block: [*c]gsl_block_int = @import("std").mem.zeroes([*c]gsl_block_int),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_int_view = extern struct {
    vector: gsl_vector_int = @import("std").mem.zeroes(gsl_vector_int),
};
pub const gsl_vector_int_view = _gsl_vector_int_view;
pub const _gsl_vector_int_const_view = extern struct {
    vector: gsl_vector_int = @import("std").mem.zeroes(gsl_vector_int),
};
pub const gsl_vector_int_const_view = _gsl_vector_int_const_view;
pub extern fn gsl_vector_int_alloc(n: usize) [*c]gsl_vector_int;
pub extern fn gsl_vector_int_calloc(n: usize) [*c]gsl_vector_int;
pub extern fn gsl_vector_int_alloc_from_block(b: [*c]gsl_block_int, offset: usize, n: usize, stride: usize) [*c]gsl_vector_int;
pub extern fn gsl_vector_int_alloc_from_vector(v: [*c]gsl_vector_int, offset: usize, n: usize, stride: usize) [*c]gsl_vector_int;
pub extern fn gsl_vector_int_free(v: [*c]gsl_vector_int) void;
pub extern fn gsl_vector_int_view_array(v: [*c]c_int, n: usize) _gsl_vector_int_view;
pub extern fn gsl_vector_int_view_array_with_stride(base: [*c]c_int, stride: usize, n: usize) _gsl_vector_int_view;
pub extern fn gsl_vector_int_const_view_array(v: [*c]const c_int, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_vector_int_const_view_array_with_stride(base: [*c]const c_int, stride: usize, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_vector_int_subvector(v: [*c]gsl_vector_int, i: usize, n: usize) _gsl_vector_int_view;
pub extern fn gsl_vector_int_subvector_with_stride(v: [*c]gsl_vector_int, i: usize, stride: usize, n: usize) _gsl_vector_int_view;
pub extern fn gsl_vector_int_const_subvector(v: [*c]const gsl_vector_int, i: usize, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_vector_int_const_subvector_with_stride(v: [*c]const gsl_vector_int, i: usize, stride: usize, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_vector_int_set_zero(v: [*c]gsl_vector_int) void;
pub extern fn gsl_vector_int_set_all(v: [*c]gsl_vector_int, x: c_int) void;
pub extern fn gsl_vector_int_set_basis(v: [*c]gsl_vector_int, i: usize) c_int;
pub extern fn gsl_vector_int_fread(stream: ?*FILE, v: [*c]gsl_vector_int) c_int;
pub extern fn gsl_vector_int_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_fscanf(stream: ?*FILE, v: [*c]gsl_vector_int) c_int;
pub extern fn gsl_vector_int_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_int, format: [*c]const u8) c_int;
pub extern fn gsl_vector_int_memcpy(dest: [*c]gsl_vector_int, src: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_reverse(v: [*c]gsl_vector_int) c_int;
pub extern fn gsl_vector_int_swap(v: [*c]gsl_vector_int, w: [*c]gsl_vector_int) c_int;
pub extern fn gsl_vector_int_swap_elements(v: [*c]gsl_vector_int, i: usize, j: usize) c_int;
pub extern fn gsl_vector_int_max(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_min(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_minmax(v: [*c]const gsl_vector_int, min_out: [*c]c_int, max_out: [*c]c_int) void;
pub extern fn gsl_vector_int_max_index(v: [*c]const gsl_vector_int) usize;
pub extern fn gsl_vector_int_min_index(v: [*c]const gsl_vector_int) usize;
pub extern fn gsl_vector_int_minmax_index(v: [*c]const gsl_vector_int, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_int_add(a: [*c]gsl_vector_int, b: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_sub(a: [*c]gsl_vector_int, b: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_mul(a: [*c]gsl_vector_int, b: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_div(a: [*c]gsl_vector_int, b: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_scale(a: [*c]gsl_vector_int, x: c_int) c_int;
pub extern fn gsl_vector_int_add_constant(a: [*c]gsl_vector_int, x: c_int) c_int;
pub extern fn gsl_vector_int_axpby(alpha: c_int, x: [*c]const gsl_vector_int, beta: c_int, y: [*c]gsl_vector_int) c_int;
pub extern fn gsl_vector_int_sum(a: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_equal(u: [*c]const gsl_vector_int, v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_isnull(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_ispos(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_isneg(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_isnonneg(v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_vector_int_get(v: [*c]const gsl_vector_int, i: usize) c_int;
pub extern fn gsl_vector_int_set(v: [*c]gsl_vector_int, i: usize, x: c_int) void;
pub extern fn gsl_vector_int_ptr(v: [*c]gsl_vector_int, i: usize) [*c]c_int;
pub extern fn gsl_vector_int_const_ptr(v: [*c]const gsl_vector_int, i: usize) [*c]const c_int;
pub const struct_gsl_block_ushort_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ushort = @import("std").mem.zeroes([*c]c_ushort),
};
pub const gsl_block_ushort = struct_gsl_block_ushort_struct;
pub extern fn gsl_block_ushort_alloc(n: usize) [*c]gsl_block_ushort;
pub extern fn gsl_block_ushort_calloc(n: usize) [*c]gsl_block_ushort;
pub extern fn gsl_block_ushort_free(b: [*c]gsl_block_ushort) void;
pub extern fn gsl_block_ushort_fread(stream: ?*FILE, b: [*c]gsl_block_ushort) c_int;
pub extern fn gsl_block_ushort_fwrite(stream: ?*FILE, b: [*c]const gsl_block_ushort) c_int;
pub extern fn gsl_block_ushort_fscanf(stream: ?*FILE, b: [*c]gsl_block_ushort) c_int;
pub extern fn gsl_block_ushort_fprintf(stream: ?*FILE, b: [*c]const gsl_block_ushort, format: [*c]const u8) c_int;
pub extern fn gsl_block_ushort_raw_fread(stream: ?*FILE, b: [*c]c_ushort, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ushort_raw_fwrite(stream: ?*FILE, b: [*c]const c_ushort, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ushort_raw_fscanf(stream: ?*FILE, b: [*c]c_ushort, n: usize, stride: usize) c_int;
pub extern fn gsl_block_ushort_raw_fprintf(stream: ?*FILE, b: [*c]const c_ushort, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_ushort_size(b: [*c]const gsl_block_ushort) usize;
pub extern fn gsl_block_ushort_data(b: [*c]const gsl_block_ushort) [*c]c_ushort;
pub const gsl_vector_ushort = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ushort = @import("std").mem.zeroes([*c]c_ushort),
    block: [*c]gsl_block_ushort = @import("std").mem.zeroes([*c]gsl_block_ushort),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_ushort_view = extern struct {
    vector: gsl_vector_ushort = @import("std").mem.zeroes(gsl_vector_ushort),
};
pub const gsl_vector_ushort_view = _gsl_vector_ushort_view;
pub const _gsl_vector_ushort_const_view = extern struct {
    vector: gsl_vector_ushort = @import("std").mem.zeroes(gsl_vector_ushort),
};
pub const gsl_vector_ushort_const_view = _gsl_vector_ushort_const_view;
pub extern fn gsl_vector_ushort_alloc(n: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_vector_ushort_calloc(n: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_vector_ushort_alloc_from_block(b: [*c]gsl_block_ushort, offset: usize, n: usize, stride: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_vector_ushort_alloc_from_vector(v: [*c]gsl_vector_ushort, offset: usize, n: usize, stride: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_vector_ushort_free(v: [*c]gsl_vector_ushort) void;
pub extern fn gsl_vector_ushort_view_array(v: [*c]c_ushort, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_vector_ushort_view_array_with_stride(base: [*c]c_ushort, stride: usize, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_vector_ushort_const_view_array(v: [*c]const c_ushort, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_vector_ushort_const_view_array_with_stride(base: [*c]const c_ushort, stride: usize, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_vector_ushort_subvector(v: [*c]gsl_vector_ushort, i: usize, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_vector_ushort_subvector_with_stride(v: [*c]gsl_vector_ushort, i: usize, stride: usize, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_vector_ushort_const_subvector(v: [*c]const gsl_vector_ushort, i: usize, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_vector_ushort_const_subvector_with_stride(v: [*c]const gsl_vector_ushort, i: usize, stride: usize, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_vector_ushort_set_zero(v: [*c]gsl_vector_ushort) void;
pub extern fn gsl_vector_ushort_set_all(v: [*c]gsl_vector_ushort, x: c_ushort) void;
pub extern fn gsl_vector_ushort_set_basis(v: [*c]gsl_vector_ushort, i: usize) c_int;
pub extern fn gsl_vector_ushort_fread(stream: ?*FILE, v: [*c]gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_fscanf(stream: ?*FILE, v: [*c]gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_ushort, format: [*c]const u8) c_int;
pub extern fn gsl_vector_ushort_memcpy(dest: [*c]gsl_vector_ushort, src: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_reverse(v: [*c]gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_swap(v: [*c]gsl_vector_ushort, w: [*c]gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_swap_elements(v: [*c]gsl_vector_ushort, i: usize, j: usize) c_int;
pub extern fn gsl_vector_ushort_max(v: [*c]const gsl_vector_ushort) c_ushort;
pub extern fn gsl_vector_ushort_min(v: [*c]const gsl_vector_ushort) c_ushort;
pub extern fn gsl_vector_ushort_minmax(v: [*c]const gsl_vector_ushort, min_out: [*c]c_ushort, max_out: [*c]c_ushort) void;
pub extern fn gsl_vector_ushort_max_index(v: [*c]const gsl_vector_ushort) usize;
pub extern fn gsl_vector_ushort_min_index(v: [*c]const gsl_vector_ushort) usize;
pub extern fn gsl_vector_ushort_minmax_index(v: [*c]const gsl_vector_ushort, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_ushort_add(a: [*c]gsl_vector_ushort, b: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_sub(a: [*c]gsl_vector_ushort, b: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_mul(a: [*c]gsl_vector_ushort, b: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_div(a: [*c]gsl_vector_ushort, b: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_scale(a: [*c]gsl_vector_ushort, x: c_ushort) c_int;
pub extern fn gsl_vector_ushort_add_constant(a: [*c]gsl_vector_ushort, x: c_ushort) c_int;
pub extern fn gsl_vector_ushort_axpby(alpha: c_ushort, x: [*c]const gsl_vector_ushort, beta: c_ushort, y: [*c]gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_sum(a: [*c]const gsl_vector_ushort) c_ushort;
pub extern fn gsl_vector_ushort_equal(u: [*c]const gsl_vector_ushort, v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_isnull(v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_ispos(v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_isneg(v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_isnonneg(v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_vector_ushort_get(v: [*c]const gsl_vector_ushort, i: usize) c_ushort;
pub extern fn gsl_vector_ushort_set(v: [*c]gsl_vector_ushort, i: usize, x: c_ushort) void;
pub extern fn gsl_vector_ushort_ptr(v: [*c]gsl_vector_ushort, i: usize) [*c]c_ushort;
pub extern fn gsl_vector_ushort_const_ptr(v: [*c]const gsl_vector_ushort, i: usize) [*c]const c_ushort;
pub const struct_gsl_block_short_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_short = @import("std").mem.zeroes([*c]c_short),
};
pub const gsl_block_short = struct_gsl_block_short_struct;
pub extern fn gsl_block_short_alloc(n: usize) [*c]gsl_block_short;
pub extern fn gsl_block_short_calloc(n: usize) [*c]gsl_block_short;
pub extern fn gsl_block_short_free(b: [*c]gsl_block_short) void;
pub extern fn gsl_block_short_fread(stream: ?*FILE, b: [*c]gsl_block_short) c_int;
pub extern fn gsl_block_short_fwrite(stream: ?*FILE, b: [*c]const gsl_block_short) c_int;
pub extern fn gsl_block_short_fscanf(stream: ?*FILE, b: [*c]gsl_block_short) c_int;
pub extern fn gsl_block_short_fprintf(stream: ?*FILE, b: [*c]const gsl_block_short, format: [*c]const u8) c_int;
pub extern fn gsl_block_short_raw_fread(stream: ?*FILE, b: [*c]c_short, n: usize, stride: usize) c_int;
pub extern fn gsl_block_short_raw_fwrite(stream: ?*FILE, b: [*c]const c_short, n: usize, stride: usize) c_int;
pub extern fn gsl_block_short_raw_fscanf(stream: ?*FILE, b: [*c]c_short, n: usize, stride: usize) c_int;
pub extern fn gsl_block_short_raw_fprintf(stream: ?*FILE, b: [*c]const c_short, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_short_size(b: [*c]const gsl_block_short) usize;
pub extern fn gsl_block_short_data(b: [*c]const gsl_block_short) [*c]c_short;
pub const gsl_vector_short = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_short = @import("std").mem.zeroes([*c]c_short),
    block: [*c]gsl_block_short = @import("std").mem.zeroes([*c]gsl_block_short),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_short_view = extern struct {
    vector: gsl_vector_short = @import("std").mem.zeroes(gsl_vector_short),
};
pub const gsl_vector_short_view = _gsl_vector_short_view;
pub const _gsl_vector_short_const_view = extern struct {
    vector: gsl_vector_short = @import("std").mem.zeroes(gsl_vector_short),
};
pub const gsl_vector_short_const_view = _gsl_vector_short_const_view;
pub extern fn gsl_vector_short_alloc(n: usize) [*c]gsl_vector_short;
pub extern fn gsl_vector_short_calloc(n: usize) [*c]gsl_vector_short;
pub extern fn gsl_vector_short_alloc_from_block(b: [*c]gsl_block_short, offset: usize, n: usize, stride: usize) [*c]gsl_vector_short;
pub extern fn gsl_vector_short_alloc_from_vector(v: [*c]gsl_vector_short, offset: usize, n: usize, stride: usize) [*c]gsl_vector_short;
pub extern fn gsl_vector_short_free(v: [*c]gsl_vector_short) void;
pub extern fn gsl_vector_short_view_array(v: [*c]c_short, n: usize) _gsl_vector_short_view;
pub extern fn gsl_vector_short_view_array_with_stride(base: [*c]c_short, stride: usize, n: usize) _gsl_vector_short_view;
pub extern fn gsl_vector_short_const_view_array(v: [*c]const c_short, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_vector_short_const_view_array_with_stride(base: [*c]const c_short, stride: usize, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_vector_short_subvector(v: [*c]gsl_vector_short, i: usize, n: usize) _gsl_vector_short_view;
pub extern fn gsl_vector_short_subvector_with_stride(v: [*c]gsl_vector_short, i: usize, stride: usize, n: usize) _gsl_vector_short_view;
pub extern fn gsl_vector_short_const_subvector(v: [*c]const gsl_vector_short, i: usize, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_vector_short_const_subvector_with_stride(v: [*c]const gsl_vector_short, i: usize, stride: usize, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_vector_short_set_zero(v: [*c]gsl_vector_short) void;
pub extern fn gsl_vector_short_set_all(v: [*c]gsl_vector_short, x: c_short) void;
pub extern fn gsl_vector_short_set_basis(v: [*c]gsl_vector_short, i: usize) c_int;
pub extern fn gsl_vector_short_fread(stream: ?*FILE, v: [*c]gsl_vector_short) c_int;
pub extern fn gsl_vector_short_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_fscanf(stream: ?*FILE, v: [*c]gsl_vector_short) c_int;
pub extern fn gsl_vector_short_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_short, format: [*c]const u8) c_int;
pub extern fn gsl_vector_short_memcpy(dest: [*c]gsl_vector_short, src: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_reverse(v: [*c]gsl_vector_short) c_int;
pub extern fn gsl_vector_short_swap(v: [*c]gsl_vector_short, w: [*c]gsl_vector_short) c_int;
pub extern fn gsl_vector_short_swap_elements(v: [*c]gsl_vector_short, i: usize, j: usize) c_int;
pub extern fn gsl_vector_short_max(v: [*c]const gsl_vector_short) c_short;
pub extern fn gsl_vector_short_min(v: [*c]const gsl_vector_short) c_short;
pub extern fn gsl_vector_short_minmax(v: [*c]const gsl_vector_short, min_out: [*c]c_short, max_out: [*c]c_short) void;
pub extern fn gsl_vector_short_max_index(v: [*c]const gsl_vector_short) usize;
pub extern fn gsl_vector_short_min_index(v: [*c]const gsl_vector_short) usize;
pub extern fn gsl_vector_short_minmax_index(v: [*c]const gsl_vector_short, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_short_add(a: [*c]gsl_vector_short, b: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_sub(a: [*c]gsl_vector_short, b: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_mul(a: [*c]gsl_vector_short, b: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_div(a: [*c]gsl_vector_short, b: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_scale(a: [*c]gsl_vector_short, x: c_short) c_int;
pub extern fn gsl_vector_short_add_constant(a: [*c]gsl_vector_short, x: c_short) c_int;
pub extern fn gsl_vector_short_axpby(alpha: c_short, x: [*c]const gsl_vector_short, beta: c_short, y: [*c]gsl_vector_short) c_int;
pub extern fn gsl_vector_short_sum(a: [*c]const gsl_vector_short) c_short;
pub extern fn gsl_vector_short_equal(u: [*c]const gsl_vector_short, v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_isnull(v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_ispos(v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_isneg(v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_isnonneg(v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_vector_short_get(v: [*c]const gsl_vector_short, i: usize) c_short;
pub extern fn gsl_vector_short_set(v: [*c]gsl_vector_short, i: usize, x: c_short) void;
pub extern fn gsl_vector_short_ptr(v: [*c]gsl_vector_short, i: usize) [*c]c_short;
pub extern fn gsl_vector_short_const_ptr(v: [*c]const gsl_vector_short, i: usize) [*c]const c_short;
pub const struct_gsl_block_uchar_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
};
pub const gsl_block_uchar = struct_gsl_block_uchar_struct;
pub extern fn gsl_block_uchar_alloc(n: usize) [*c]gsl_block_uchar;
pub extern fn gsl_block_uchar_calloc(n: usize) [*c]gsl_block_uchar;
pub extern fn gsl_block_uchar_free(b: [*c]gsl_block_uchar) void;
pub extern fn gsl_block_uchar_fread(stream: ?*FILE, b: [*c]gsl_block_uchar) c_int;
pub extern fn gsl_block_uchar_fwrite(stream: ?*FILE, b: [*c]const gsl_block_uchar) c_int;
pub extern fn gsl_block_uchar_fscanf(stream: ?*FILE, b: [*c]gsl_block_uchar) c_int;
pub extern fn gsl_block_uchar_fprintf(stream: ?*FILE, b: [*c]const gsl_block_uchar, format: [*c]const u8) c_int;
pub extern fn gsl_block_uchar_raw_fread(stream: ?*FILE, b: [*c]u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uchar_raw_fwrite(stream: ?*FILE, b: [*c]const u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uchar_raw_fscanf(stream: ?*FILE, b: [*c]u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_uchar_raw_fprintf(stream: ?*FILE, b: [*c]const u8, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_uchar_size(b: [*c]const gsl_block_uchar) usize;
pub extern fn gsl_block_uchar_data(b: [*c]const gsl_block_uchar) [*c]u8;
pub const gsl_vector_uchar = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
    block: [*c]gsl_block_uchar = @import("std").mem.zeroes([*c]gsl_block_uchar),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_uchar_view = extern struct {
    vector: gsl_vector_uchar = @import("std").mem.zeroes(gsl_vector_uchar),
};
pub const gsl_vector_uchar_view = _gsl_vector_uchar_view;
pub const _gsl_vector_uchar_const_view = extern struct {
    vector: gsl_vector_uchar = @import("std").mem.zeroes(gsl_vector_uchar),
};
pub const gsl_vector_uchar_const_view = _gsl_vector_uchar_const_view;
pub extern fn gsl_vector_uchar_alloc(n: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_vector_uchar_calloc(n: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_vector_uchar_alloc_from_block(b: [*c]gsl_block_uchar, offset: usize, n: usize, stride: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_vector_uchar_alloc_from_vector(v: [*c]gsl_vector_uchar, offset: usize, n: usize, stride: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_vector_uchar_free(v: [*c]gsl_vector_uchar) void;
pub extern fn gsl_vector_uchar_view_array(v: [*c]u8, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_vector_uchar_view_array_with_stride(base: [*c]u8, stride: usize, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_vector_uchar_const_view_array(v: [*c]const u8, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_vector_uchar_const_view_array_with_stride(base: [*c]const u8, stride: usize, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_vector_uchar_subvector(v: [*c]gsl_vector_uchar, i: usize, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_vector_uchar_subvector_with_stride(v: [*c]gsl_vector_uchar, i: usize, stride: usize, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_vector_uchar_const_subvector(v: [*c]const gsl_vector_uchar, i: usize, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_vector_uchar_const_subvector_with_stride(v: [*c]const gsl_vector_uchar, i: usize, stride: usize, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_vector_uchar_set_zero(v: [*c]gsl_vector_uchar) void;
pub extern fn gsl_vector_uchar_set_all(v: [*c]gsl_vector_uchar, x: u8) void;
pub extern fn gsl_vector_uchar_set_basis(v: [*c]gsl_vector_uchar, i: usize) c_int;
pub extern fn gsl_vector_uchar_fread(stream: ?*FILE, v: [*c]gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_fscanf(stream: ?*FILE, v: [*c]gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_uchar, format: [*c]const u8) c_int;
pub extern fn gsl_vector_uchar_memcpy(dest: [*c]gsl_vector_uchar, src: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_reverse(v: [*c]gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_swap(v: [*c]gsl_vector_uchar, w: [*c]gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_swap_elements(v: [*c]gsl_vector_uchar, i: usize, j: usize) c_int;
pub extern fn gsl_vector_uchar_max(v: [*c]const gsl_vector_uchar) u8;
pub extern fn gsl_vector_uchar_min(v: [*c]const gsl_vector_uchar) u8;
pub extern fn gsl_vector_uchar_minmax(v: [*c]const gsl_vector_uchar, min_out: [*c]u8, max_out: [*c]u8) void;
pub extern fn gsl_vector_uchar_max_index(v: [*c]const gsl_vector_uchar) usize;
pub extern fn gsl_vector_uchar_min_index(v: [*c]const gsl_vector_uchar) usize;
pub extern fn gsl_vector_uchar_minmax_index(v: [*c]const gsl_vector_uchar, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_uchar_add(a: [*c]gsl_vector_uchar, b: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_sub(a: [*c]gsl_vector_uchar, b: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_mul(a: [*c]gsl_vector_uchar, b: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_div(a: [*c]gsl_vector_uchar, b: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_scale(a: [*c]gsl_vector_uchar, x: u8) c_int;
pub extern fn gsl_vector_uchar_add_constant(a: [*c]gsl_vector_uchar, x: u8) c_int;
pub extern fn gsl_vector_uchar_axpby(alpha: u8, x: [*c]const gsl_vector_uchar, beta: u8, y: [*c]gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_sum(a: [*c]const gsl_vector_uchar) u8;
pub extern fn gsl_vector_uchar_equal(u: [*c]const gsl_vector_uchar, v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_isnull(v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_ispos(v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_isneg(v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_isnonneg(v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_vector_uchar_get(v: [*c]const gsl_vector_uchar, i: usize) u8;
pub extern fn gsl_vector_uchar_set(v: [*c]gsl_vector_uchar, i: usize, x: u8) void;
pub extern fn gsl_vector_uchar_ptr(v: [*c]gsl_vector_uchar, i: usize) [*c]u8;
pub extern fn gsl_vector_uchar_const_ptr(v: [*c]const gsl_vector_uchar, i: usize) [*c]const u8;
pub const struct_gsl_block_char_struct = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
};
pub const gsl_block_char = struct_gsl_block_char_struct;
pub extern fn gsl_block_char_alloc(n: usize) [*c]gsl_block_char;
pub extern fn gsl_block_char_calloc(n: usize) [*c]gsl_block_char;
pub extern fn gsl_block_char_free(b: [*c]gsl_block_char) void;
pub extern fn gsl_block_char_fread(stream: ?*FILE, b: [*c]gsl_block_char) c_int;
pub extern fn gsl_block_char_fwrite(stream: ?*FILE, b: [*c]const gsl_block_char) c_int;
pub extern fn gsl_block_char_fscanf(stream: ?*FILE, b: [*c]gsl_block_char) c_int;
pub extern fn gsl_block_char_fprintf(stream: ?*FILE, b: [*c]const gsl_block_char, format: [*c]const u8) c_int;
pub extern fn gsl_block_char_raw_fread(stream: ?*FILE, b: [*c]u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_char_raw_fwrite(stream: ?*FILE, b: [*c]const u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_char_raw_fscanf(stream: ?*FILE, b: [*c]u8, n: usize, stride: usize) c_int;
pub extern fn gsl_block_char_raw_fprintf(stream: ?*FILE, b: [*c]const u8, n: usize, stride: usize, format: [*c]const u8) c_int;
pub extern fn gsl_block_char_size(b: [*c]const gsl_block_char) usize;
pub extern fn gsl_block_char_data(b: [*c]const gsl_block_char) [*c]u8;
pub const gsl_vector_char = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    stride: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
    block: [*c]gsl_block_char = @import("std").mem.zeroes([*c]gsl_block_char),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_vector_char_view = extern struct {
    vector: gsl_vector_char = @import("std").mem.zeroes(gsl_vector_char),
};
pub const gsl_vector_char_view = _gsl_vector_char_view;
pub const _gsl_vector_char_const_view = extern struct {
    vector: gsl_vector_char = @import("std").mem.zeroes(gsl_vector_char),
};
pub const gsl_vector_char_const_view = _gsl_vector_char_const_view;
pub extern fn gsl_vector_char_alloc(n: usize) [*c]gsl_vector_char;
pub extern fn gsl_vector_char_calloc(n: usize) [*c]gsl_vector_char;
pub extern fn gsl_vector_char_alloc_from_block(b: [*c]gsl_block_char, offset: usize, n: usize, stride: usize) [*c]gsl_vector_char;
pub extern fn gsl_vector_char_alloc_from_vector(v: [*c]gsl_vector_char, offset: usize, n: usize, stride: usize) [*c]gsl_vector_char;
pub extern fn gsl_vector_char_free(v: [*c]gsl_vector_char) void;
pub extern fn gsl_vector_char_view_array(v: [*c]u8, n: usize) _gsl_vector_char_view;
pub extern fn gsl_vector_char_view_array_with_stride(base: [*c]u8, stride: usize, n: usize) _gsl_vector_char_view;
pub extern fn gsl_vector_char_const_view_array(v: [*c]const u8, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_vector_char_const_view_array_with_stride(base: [*c]const u8, stride: usize, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_vector_char_subvector(v: [*c]gsl_vector_char, i: usize, n: usize) _gsl_vector_char_view;
pub extern fn gsl_vector_char_subvector_with_stride(v: [*c]gsl_vector_char, i: usize, stride: usize, n: usize) _gsl_vector_char_view;
pub extern fn gsl_vector_char_const_subvector(v: [*c]const gsl_vector_char, i: usize, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_vector_char_const_subvector_with_stride(v: [*c]const gsl_vector_char, i: usize, stride: usize, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_vector_char_set_zero(v: [*c]gsl_vector_char) void;
pub extern fn gsl_vector_char_set_all(v: [*c]gsl_vector_char, x: u8) void;
pub extern fn gsl_vector_char_set_basis(v: [*c]gsl_vector_char, i: usize) c_int;
pub extern fn gsl_vector_char_fread(stream: ?*FILE, v: [*c]gsl_vector_char) c_int;
pub extern fn gsl_vector_char_fwrite(stream: ?*FILE, v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_fscanf(stream: ?*FILE, v: [*c]gsl_vector_char) c_int;
pub extern fn gsl_vector_char_fprintf(stream: ?*FILE, v: [*c]const gsl_vector_char, format: [*c]const u8) c_int;
pub extern fn gsl_vector_char_memcpy(dest: [*c]gsl_vector_char, src: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_reverse(v: [*c]gsl_vector_char) c_int;
pub extern fn gsl_vector_char_swap(v: [*c]gsl_vector_char, w: [*c]gsl_vector_char) c_int;
pub extern fn gsl_vector_char_swap_elements(v: [*c]gsl_vector_char, i: usize, j: usize) c_int;
pub extern fn gsl_vector_char_max(v: [*c]const gsl_vector_char) u8;
pub extern fn gsl_vector_char_min(v: [*c]const gsl_vector_char) u8;
pub extern fn gsl_vector_char_minmax(v: [*c]const gsl_vector_char, min_out: [*c]u8, max_out: [*c]u8) void;
pub extern fn gsl_vector_char_max_index(v: [*c]const gsl_vector_char) usize;
pub extern fn gsl_vector_char_min_index(v: [*c]const gsl_vector_char) usize;
pub extern fn gsl_vector_char_minmax_index(v: [*c]const gsl_vector_char, imin: [*c]usize, imax: [*c]usize) void;
pub extern fn gsl_vector_char_add(a: [*c]gsl_vector_char, b: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_sub(a: [*c]gsl_vector_char, b: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_mul(a: [*c]gsl_vector_char, b: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_div(a: [*c]gsl_vector_char, b: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_scale(a: [*c]gsl_vector_char, x: u8) c_int;
pub extern fn gsl_vector_char_add_constant(a: [*c]gsl_vector_char, x: u8) c_int;
pub extern fn gsl_vector_char_axpby(alpha: u8, x: [*c]const gsl_vector_char, beta: u8, y: [*c]gsl_vector_char) c_int;
pub extern fn gsl_vector_char_sum(a: [*c]const gsl_vector_char) u8;
pub extern fn gsl_vector_char_equal(u: [*c]const gsl_vector_char, v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_isnull(v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_ispos(v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_isneg(v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_isnonneg(v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_vector_char_get(v: [*c]const gsl_vector_char, i: usize) u8;
pub extern fn gsl_vector_char_set(v: [*c]gsl_vector_char, i: usize, x: u8) void;
pub extern fn gsl_vector_char_ptr(v: [*c]gsl_vector_char, i: usize) [*c]u8;
pub extern fn gsl_vector_char_const_ptr(v: [*c]const gsl_vector_char, i: usize) [*c]const u8;
pub const ptrdiff_t = c_long;
pub const max_align_t = extern struct {
    __clang_max_align_nonce1: c_longlong align(8) = @import("std").mem.zeroes(c_longlong),
    __clang_max_align_nonce2: c_longdouble align(16) = @import("std").mem.zeroes(c_longdouble),
};
pub const CblasRowMajor: c_int = 101;
pub const CblasColMajor: c_int = 102;
pub const enum_CBLAS_ORDER = c_uint;
pub const CblasNoTrans: c_int = 111;
pub const CblasTrans: c_int = 112;
pub const CblasConjTrans: c_int = 113;
pub const enum_CBLAS_TRANSPOSE = c_uint;
pub const CblasUpper: c_int = 121;
pub const CblasLower: c_int = 122;
pub const enum_CBLAS_UPLO = c_uint;
pub const CblasNonUnit: c_int = 131;
pub const CblasUnit: c_int = 132;
pub const enum_CBLAS_DIAG = c_uint;
pub const CblasLeft: c_int = 141;
pub const CblasRight: c_int = 142;
pub const enum_CBLAS_SIDE = c_uint;
pub extern fn cblas_sdsdot(N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int) f32;
pub extern fn cblas_dsdot(N: c_int, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int) f64;
pub extern fn cblas_sdot(N: c_int, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int) f32;
pub extern fn cblas_ddot(N: c_int, X: [*c]const f64, incX: c_int, Y: [*c]const f64, incY: c_int) f64;
pub extern fn cblas_cdotu_sub(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, dotu: ?*anyopaque) void;
pub extern fn cblas_cdotc_sub(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, dotc: ?*anyopaque) void;
pub extern fn cblas_zdotu_sub(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, dotu: ?*anyopaque) void;
pub extern fn cblas_zdotc_sub(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, dotc: ?*anyopaque) void;
pub extern fn cblas_snrm2(N: c_int, X: [*c]const f32, incX: c_int) f32;
pub extern fn cblas_sasum(N: c_int, X: [*c]const f32, incX: c_int) f32;
pub extern fn cblas_dnrm2(N: c_int, X: [*c]const f64, incX: c_int) f64;
pub extern fn cblas_dasum(N: c_int, X: [*c]const f64, incX: c_int) f64;
pub extern fn cblas_scnrm2(N: c_int, X: ?*const anyopaque, incX: c_int) f32;
pub extern fn cblas_scasum(N: c_int, X: ?*const anyopaque, incX: c_int) f32;
pub extern fn cblas_dznrm2(N: c_int, X: ?*const anyopaque, incX: c_int) f64;
pub extern fn cblas_dzasum(N: c_int, X: ?*const anyopaque, incX: c_int) f64;
pub extern fn cblas_isamax(N: c_int, X: [*c]const f32, incX: c_int) usize;
pub extern fn cblas_idamax(N: c_int, X: [*c]const f64, incX: c_int) usize;
pub extern fn cblas_icamax(N: c_int, X: ?*const anyopaque, incX: c_int) usize;
pub extern fn cblas_izamax(N: c_int, X: ?*const anyopaque, incX: c_int) usize;
pub extern fn cblas_sswap(N: c_int, X: [*c]f32, incX: c_int, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_scopy(N: c_int, X: [*c]const f32, incX: c_int, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_saxpy(N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_dswap(N: c_int, X: [*c]f64, incX: c_int, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dcopy(N: c_int, X: [*c]const f64, incX: c_int, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_daxpy(N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_cswap(N: c_int, X: ?*anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_ccopy(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_caxpy(N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zswap(N: c_int, X: ?*anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zcopy(N: c_int, X: ?*const anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zaxpy(N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_srotg(a: [*c]f32, b: [*c]f32, c: [*c]f32, s: [*c]f32) void;
pub extern fn cblas_srotmg(d1: [*c]f32, d2: [*c]f32, b1: [*c]f32, b2: f32, P: [*c]f32) void;
pub extern fn cblas_srot(N: c_int, X: [*c]f32, incX: c_int, Y: [*c]f32, incY: c_int, c: f32, s: f32) void;
pub extern fn cblas_srotm(N: c_int, X: [*c]f32, incX: c_int, Y: [*c]f32, incY: c_int, P: [*c]const f32) void;
pub extern fn cblas_drotg(a: [*c]f64, b: [*c]f64, c: [*c]f64, s: [*c]f64) void;
pub extern fn cblas_drotmg(d1: [*c]f64, d2: [*c]f64, b1: [*c]f64, b2: f64, P: [*c]f64) void;
pub extern fn cblas_drot(N: c_int, X: [*c]f64, incX: c_int, Y: [*c]f64, incY: c_int, c: f64, s: f64) void;
pub extern fn cblas_drotm(N: c_int, X: [*c]f64, incX: c_int, Y: [*c]f64, incY: c_int, P: [*c]const f64) void;
pub extern fn cblas_sscal(N: c_int, alpha: f32, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_dscal(N: c_int, alpha: f64, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_cscal(N: c_int, alpha: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_zscal(N: c_int, alpha: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_csscal(N: c_int, alpha: f32, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_zdscal(N: c_int, alpha: f64, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_sgemv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, alpha: f32, A: [*c]const f32, lda: c_int, X: [*c]const f32, incX: c_int, beta: f32, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_sgbmv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, KL: c_int, KU: c_int, alpha: f32, A: [*c]const f32, lda: c_int, X: [*c]const f32, incX: c_int, beta: f32, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_strmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: [*c]const f32, lda: c_int, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_stbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: [*c]const f32, lda: c_int, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_stpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: [*c]const f32, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_strsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: [*c]const f32, lda: c_int, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_stbsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: [*c]const f32, lda: c_int, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_stpsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: [*c]const f32, X: [*c]f32, incX: c_int) void;
pub extern fn cblas_dgemv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, alpha: f64, A: [*c]const f64, lda: c_int, X: [*c]const f64, incX: c_int, beta: f64, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dgbmv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, KL: c_int, KU: c_int, alpha: f64, A: [*c]const f64, lda: c_int, X: [*c]const f64, incX: c_int, beta: f64, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dtrmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: [*c]const f64, lda: c_int, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_dtbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: [*c]const f64, lda: c_int, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_dtpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: [*c]const f64, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_dtrsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: [*c]const f64, lda: c_int, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_dtbsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: [*c]const f64, lda: c_int, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_dtpsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: [*c]const f64, X: [*c]f64, incX: c_int) void;
pub extern fn cblas_cgemv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_cgbmv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, KL: c_int, KU: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_ctrmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ctbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ctpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ctrsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ctbsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ctpsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_zgemv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zgbmv(order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, KL: c_int, KU: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_ztrmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ztbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ztpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ztrsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ztbsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, K: c_int, A: ?*const anyopaque, lda: c_int, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ztpsv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, N: c_int, Ap: ?*const anyopaque, X: ?*anyopaque, incX: c_int) void;
pub extern fn cblas_ssymv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, A: [*c]const f32, lda: c_int, X: [*c]const f32, incX: c_int, beta: f32, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_ssbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, K: c_int, alpha: f32, A: [*c]const f32, lda: c_int, X: [*c]const f32, incX: c_int, beta: f32, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_sspmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, Ap: [*c]const f32, X: [*c]const f32, incX: c_int, beta: f32, Y: [*c]f32, incY: c_int) void;
pub extern fn cblas_sger(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int, A: [*c]f32, lda: c_int) void;
pub extern fn cblas_ssyr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, A: [*c]f32, lda: c_int) void;
pub extern fn cblas_sspr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Ap: [*c]f32) void;
pub extern fn cblas_ssyr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int, A: [*c]f32, lda: c_int) void;
pub extern fn cblas_sspr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: [*c]const f32, incX: c_int, Y: [*c]const f32, incY: c_int, A: [*c]f32) void;
pub extern fn cblas_dsymv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, A: [*c]const f64, lda: c_int, X: [*c]const f64, incX: c_int, beta: f64, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dsbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, K: c_int, alpha: f64, A: [*c]const f64, lda: c_int, X: [*c]const f64, incX: c_int, beta: f64, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dspmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, Ap: [*c]const f64, X: [*c]const f64, incX: c_int, beta: f64, Y: [*c]f64, incY: c_int) void;
pub extern fn cblas_dger(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, Y: [*c]const f64, incY: c_int, A: [*c]f64, lda: c_int) void;
pub extern fn cblas_dsyr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, A: [*c]f64, lda: c_int) void;
pub extern fn cblas_dspr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, Ap: [*c]f64) void;
pub extern fn cblas_dsyr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, Y: [*c]const f64, incY: c_int, A: [*c]f64, lda: c_int) void;
pub extern fn cblas_dspr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: [*c]const f64, incX: c_int, Y: [*c]const f64, incY: c_int, A: [*c]f64) void;
pub extern fn cblas_chemv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_chbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_chpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, Ap: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_cgeru(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_cgerc(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_cher(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: ?*const anyopaque, incX: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_chpr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f32, X: ?*const anyopaque, incX: c_int, A: ?*anyopaque) void;
pub extern fn cblas_cher2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_chpr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, Ap: ?*anyopaque) void;
pub extern fn cblas_zhemv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zhbmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zhpmv(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, Ap: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, beta: ?*const anyopaque, Y: ?*anyopaque, incY: c_int) void;
pub extern fn cblas_zgeru(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_zgerc(order: enum_CBLAS_ORDER, M: c_int, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_zher(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: ?*const anyopaque, incX: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_zhpr(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: f64, X: ?*const anyopaque, incX: c_int, A: ?*anyopaque) void;
pub extern fn cblas_zher2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, A: ?*anyopaque, lda: c_int) void;
pub extern fn cblas_zhpr2(order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, N: c_int, alpha: ?*const anyopaque, X: ?*const anyopaque, incX: c_int, Y: ?*const anyopaque, incY: c_int, Ap: ?*anyopaque) void;
pub extern fn cblas_sgemm(Order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, TransB: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, K: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]const f32, ldb: c_int, beta: f32, C: [*c]f32, ldc: c_int) void;
pub extern fn cblas_ssymm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]const f32, ldb: c_int, beta: f32, C: [*c]f32, ldc: c_int) void;
pub extern fn cblas_ssyrk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f32, A: [*c]const f32, lda: c_int, beta: f32, C: [*c]f32, ldc: c_int) void;
pub extern fn cblas_ssyr2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]const f32, ldb: c_int, beta: f32, C: [*c]f32, ldc: c_int) void;
pub extern fn cblas_strmm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]f32, ldb: c_int) void;
pub extern fn cblas_strsm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]f32, ldb: c_int) void;
pub extern fn cblas_dgemm(Order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, TransB: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, K: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]const f64, ldb: c_int, beta: f64, C: [*c]f64, ldc: c_int) void;
pub extern fn cblas_dsymm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]const f64, ldb: c_int, beta: f64, C: [*c]f64, ldc: c_int) void;
pub extern fn cblas_dsyrk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f64, A: [*c]const f64, lda: c_int, beta: f64, C: [*c]f64, ldc: c_int) void;
pub extern fn cblas_dsyr2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]const f64, ldb: c_int, beta: f64, C: [*c]f64, ldc: c_int) void;
pub extern fn cblas_dtrmm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]f64, ldb: c_int) void;
pub extern fn cblas_dtrsm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]f64, ldb: c_int) void;
pub extern fn cblas_cgemm(Order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, TransB: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_csymm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_csyrk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_csyr2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_ctrmm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*anyopaque, ldb: c_int) void;
pub extern fn cblas_ctrsm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*anyopaque, ldb: c_int) void;
pub extern fn cblas_zgemm(Order: enum_CBLAS_ORDER, TransA: enum_CBLAS_TRANSPOSE, TransB: enum_CBLAS_TRANSPOSE, M: c_int, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zsymm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zsyrk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zsyr2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_ztrmm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*anyopaque, ldb: c_int) void;
pub extern fn cblas_ztrsm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, TransA: enum_CBLAS_TRANSPOSE, Diag: enum_CBLAS_DIAG, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*anyopaque, ldb: c_int) void;
pub extern fn cblas_chemm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_cherk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f32, A: ?*const anyopaque, lda: c_int, beta: f32, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_cher2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: f32, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zhemm(Order: enum_CBLAS_ORDER, Side: enum_CBLAS_SIDE, Uplo: enum_CBLAS_UPLO, M: c_int, N: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: ?*const anyopaque, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zherk(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: f64, A: ?*const anyopaque, lda: c_int, beta: f64, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_zher2k(Order: enum_CBLAS_ORDER, Uplo: enum_CBLAS_UPLO, Trans: enum_CBLAS_TRANSPOSE, N: c_int, K: c_int, alpha: ?*const anyopaque, A: ?*const anyopaque, lda: c_int, B: ?*const anyopaque, ldb: c_int, beta: f64, C: ?*anyopaque, ldc: c_int) void;
pub extern fn cblas_xerbla(p: c_int, rout: [*c]const u8, form: [*c]const u8, ...) void;
pub const CBLAS_INDEX_t = usize;
pub const CBLAS_ORDER_t = enum_CBLAS_ORDER;
pub const CBLAS_TRANSPOSE_t = enum_CBLAS_TRANSPOSE;
pub const CBLAS_UPLO_t = enum_CBLAS_UPLO;
pub const CBLAS_DIAG_t = enum_CBLAS_DIAG;
pub const CBLAS_SIDE_t = enum_CBLAS_SIDE;
pub const gsl_matrix_complex_long_double = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
    block: [*c]gsl_block_complex_long_double = @import("std").mem.zeroes([*c]gsl_block_complex_long_double),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_complex_long_double_view = extern struct {
    matrix: gsl_matrix_complex_long_double = @import("std").mem.zeroes(gsl_matrix_complex_long_double),
};
pub const gsl_matrix_complex_long_double_view = _gsl_matrix_complex_long_double_view;
pub const _gsl_matrix_complex_long_double_const_view = extern struct {
    matrix: gsl_matrix_complex_long_double = @import("std").mem.zeroes(gsl_matrix_complex_long_double),
};
pub const gsl_matrix_complex_long_double_const_view = _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_alloc(n1: usize, n2: usize) [*c]gsl_matrix_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_calloc(n1: usize, n2: usize) [*c]gsl_matrix_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_alloc_from_block(b: [*c]gsl_block_complex_long_double, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_alloc_from_matrix(b: [*c]gsl_matrix_complex_long_double, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_complex_long_double;
pub extern fn gsl_vector_complex_long_double_alloc_row_from_matrix(m: [*c]gsl_matrix_complex_long_double, i: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_vector_complex_long_double_alloc_col_from_matrix(m: [*c]gsl_matrix_complex_long_double, j: usize) [*c]gsl_vector_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_free(m: [*c]gsl_matrix_complex_long_double) void;
pub extern fn gsl_matrix_complex_long_double_submatrix(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_row(m: [*c]gsl_matrix_complex_long_double, i: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_column(m: [*c]gsl_matrix_complex_long_double, j: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_diagonal(m: [*c]gsl_matrix_complex_long_double) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_subdiagonal(m: [*c]gsl_matrix_complex_long_double, k: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_superdiagonal(m: [*c]gsl_matrix_complex_long_double, k: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_subrow(m: [*c]gsl_matrix_complex_long_double, i: usize, offset: usize, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_subcolumn(m: [*c]gsl_matrix_complex_long_double, j: usize, offset: usize, n: usize) _gsl_vector_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_view_array(base: [*c]c_longdouble, n1: usize, n2: usize) _gsl_matrix_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_view_array_with_tda(base: [*c]c_longdouble, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_view_vector(v: [*c]gsl_vector_complex_long_double, n1: usize, n2: usize) _gsl_matrix_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_view_vector_with_tda(v: [*c]gsl_vector_complex_long_double, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_long_double_view;
pub extern fn gsl_matrix_complex_long_double_const_submatrix(m: [*c]const gsl_matrix_complex_long_double, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_row(m: [*c]const gsl_matrix_complex_long_double, i: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_column(m: [*c]const gsl_matrix_complex_long_double, j: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_diagonal(m: [*c]const gsl_matrix_complex_long_double) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_subdiagonal(m: [*c]const gsl_matrix_complex_long_double, k: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_superdiagonal(m: [*c]const gsl_matrix_complex_long_double, k: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_subrow(m: [*c]const gsl_matrix_complex_long_double, i: usize, offset: usize, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_subcolumn(m: [*c]const gsl_matrix_complex_long_double, j: usize, offset: usize, n: usize) _gsl_vector_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_view_array(base: [*c]const c_longdouble, n1: usize, n2: usize) _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_view_array_with_tda(base: [*c]const c_longdouble, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_view_vector(v: [*c]const gsl_vector_complex_long_double, n1: usize, n2: usize) _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_const_view_vector_with_tda(v: [*c]const gsl_vector_complex_long_double, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_long_double_const_view;
pub extern fn gsl_matrix_complex_long_double_set_zero(m: [*c]gsl_matrix_complex_long_double) void;
pub extern fn gsl_matrix_complex_long_double_set_identity(m: [*c]gsl_matrix_complex_long_double) void;
pub extern fn gsl_matrix_complex_long_double_set_all(m: [*c]gsl_matrix_complex_long_double, x: gsl_complex_long_double) void;
pub extern fn gsl_matrix_complex_long_double_fread(stream: ?*FILE, m: [*c]gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_complex_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_complex_long_double_memcpy(dest: [*c]gsl_matrix_complex_long_double, src: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_swap(m1: [*c]gsl_matrix_complex_long_double, m2: [*c]gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex_long_double, src: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_swap_rows(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_long_double_swap_columns(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_long_double_swap_rowcol(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_long_double_transpose(m: [*c]gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_transpose_memcpy(dest: [*c]gsl_matrix_complex_long_double, src: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex_long_double, src: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_conjtrans_memcpy(dest: [*c]gsl_matrix_complex_long_double, src: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_equal(a: [*c]const gsl_matrix_complex_long_double, b: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_isnull(m: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_ispos(m: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_isneg(m: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_isnonneg(m: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_add(a: [*c]gsl_matrix_complex_long_double, b: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_sub(a: [*c]gsl_matrix_complex_long_double, b: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_mul_elements(a: [*c]gsl_matrix_complex_long_double, b: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_div_elements(a: [*c]gsl_matrix_complex_long_double, b: [*c]const gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_scale(a: [*c]gsl_matrix_complex_long_double, x: gsl_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_scale_rows(a: [*c]gsl_matrix_complex_long_double, x: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_scale_columns(a: [*c]gsl_matrix_complex_long_double, x: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_add_constant(a: [*c]gsl_matrix_complex_long_double, x: gsl_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_add_diagonal(a: [*c]gsl_matrix_complex_long_double, x: gsl_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_conjugate(a: [*c]gsl_matrix_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_get_row(v: [*c]gsl_vector_complex_long_double, m: [*c]const gsl_matrix_complex_long_double, i: usize) c_int;
pub extern fn gsl_matrix_complex_long_double_get_col(v: [*c]gsl_vector_complex_long_double, m: [*c]const gsl_matrix_complex_long_double, j: usize) c_int;
pub extern fn gsl_matrix_complex_long_double_set_row(m: [*c]gsl_matrix_complex_long_double, i: usize, v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_set_col(m: [*c]gsl_matrix_complex_long_double, j: usize, v: [*c]const gsl_vector_complex_long_double) c_int;
pub extern fn gsl_matrix_complex_long_double_get(m: [*c]const gsl_matrix_complex_long_double, i: usize, j: usize) gsl_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_set(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize, x: gsl_complex_long_double) void;
pub extern fn gsl_matrix_complex_long_double_ptr(m: [*c]gsl_matrix_complex_long_double, i: usize, j: usize) [*c]gsl_complex_long_double;
pub extern fn gsl_matrix_complex_long_double_const_ptr(m: [*c]const gsl_matrix_complex_long_double, i: usize, j: usize) [*c]const gsl_complex_long_double;
pub const gsl_matrix_complex = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    block: [*c]gsl_block_complex = @import("std").mem.zeroes([*c]gsl_block_complex),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_complex_view = extern struct {
    matrix: gsl_matrix_complex = @import("std").mem.zeroes(gsl_matrix_complex),
};
pub const gsl_matrix_complex_view = _gsl_matrix_complex_view;
pub const _gsl_matrix_complex_const_view = extern struct {
    matrix: gsl_matrix_complex = @import("std").mem.zeroes(gsl_matrix_complex),
};
pub const gsl_matrix_complex_const_view = _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_alloc(n1: usize, n2: usize) [*c]gsl_matrix_complex;
pub extern fn gsl_matrix_complex_calloc(n1: usize, n2: usize) [*c]gsl_matrix_complex;
pub extern fn gsl_matrix_complex_alloc_from_block(b: [*c]gsl_block_complex, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_complex;
pub extern fn gsl_matrix_complex_alloc_from_matrix(b: [*c]gsl_matrix_complex, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_complex;
pub extern fn gsl_vector_complex_alloc_row_from_matrix(m: [*c]gsl_matrix_complex, i: usize) [*c]gsl_vector_complex;
pub extern fn gsl_vector_complex_alloc_col_from_matrix(m: [*c]gsl_matrix_complex, j: usize) [*c]gsl_vector_complex;
pub extern fn gsl_matrix_complex_free(m: [*c]gsl_matrix_complex) void;
pub extern fn gsl_matrix_complex_submatrix(m: [*c]gsl_matrix_complex, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_view;
pub extern fn gsl_matrix_complex_row(m: [*c]gsl_matrix_complex, i: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_column(m: [*c]gsl_matrix_complex, j: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_diagonal(m: [*c]gsl_matrix_complex) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_subdiagonal(m: [*c]gsl_matrix_complex, k: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_superdiagonal(m: [*c]gsl_matrix_complex, k: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_subrow(m: [*c]gsl_matrix_complex, i: usize, offset: usize, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_subcolumn(m: [*c]gsl_matrix_complex, j: usize, offset: usize, n: usize) _gsl_vector_complex_view;
pub extern fn gsl_matrix_complex_view_array(base: [*c]f64, n1: usize, n2: usize) _gsl_matrix_complex_view;
pub extern fn gsl_matrix_complex_view_array_with_tda(base: [*c]f64, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_view;
pub extern fn gsl_matrix_complex_view_vector(v: [*c]gsl_vector_complex, n1: usize, n2: usize) _gsl_matrix_complex_view;
pub extern fn gsl_matrix_complex_view_vector_with_tda(v: [*c]gsl_vector_complex, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_view;
pub extern fn gsl_matrix_complex_const_submatrix(m: [*c]const gsl_matrix_complex, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_const_row(m: [*c]const gsl_matrix_complex, i: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_column(m: [*c]const gsl_matrix_complex, j: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_diagonal(m: [*c]const gsl_matrix_complex) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_subdiagonal(m: [*c]const gsl_matrix_complex, k: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_superdiagonal(m: [*c]const gsl_matrix_complex, k: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_subrow(m: [*c]const gsl_matrix_complex, i: usize, offset: usize, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_subcolumn(m: [*c]const gsl_matrix_complex, j: usize, offset: usize, n: usize) _gsl_vector_complex_const_view;
pub extern fn gsl_matrix_complex_const_view_array(base: [*c]const f64, n1: usize, n2: usize) _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_const_view_array_with_tda(base: [*c]const f64, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_const_view_vector(v: [*c]const gsl_vector_complex, n1: usize, n2: usize) _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_const_view_vector_with_tda(v: [*c]const gsl_vector_complex, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_const_view;
pub extern fn gsl_matrix_complex_set_zero(m: [*c]gsl_matrix_complex) void;
pub extern fn gsl_matrix_complex_set_identity(m: [*c]gsl_matrix_complex) void;
pub extern fn gsl_matrix_complex_set_all(m: [*c]gsl_matrix_complex, x: gsl_complex) void;
pub extern fn gsl_matrix_complex_fread(stream: ?*FILE, m: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_complex, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_complex_memcpy(dest: [*c]gsl_matrix_complex, src: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_swap(m1: [*c]gsl_matrix_complex, m2: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex, src: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_swap_rows(m: [*c]gsl_matrix_complex, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_swap_columns(m: [*c]gsl_matrix_complex, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_swap_rowcol(m: [*c]gsl_matrix_complex, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_transpose(m: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_transpose_memcpy(dest: [*c]gsl_matrix_complex, src: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex, src: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_conjtrans_memcpy(dest: [*c]gsl_matrix_complex, src: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_equal(a: [*c]const gsl_matrix_complex, b: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_isnull(m: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_ispos(m: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_isneg(m: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_isnonneg(m: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_add(a: [*c]gsl_matrix_complex, b: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_sub(a: [*c]gsl_matrix_complex, b: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_mul_elements(a: [*c]gsl_matrix_complex, b: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_div_elements(a: [*c]gsl_matrix_complex, b: [*c]const gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_scale(a: [*c]gsl_matrix_complex, x: gsl_complex) c_int;
pub extern fn gsl_matrix_complex_scale_rows(a: [*c]gsl_matrix_complex, x: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_matrix_complex_scale_columns(a: [*c]gsl_matrix_complex, x: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_matrix_complex_add_constant(a: [*c]gsl_matrix_complex, x: gsl_complex) c_int;
pub extern fn gsl_matrix_complex_add_diagonal(a: [*c]gsl_matrix_complex, x: gsl_complex) c_int;
pub extern fn gsl_matrix_complex_conjugate(a: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_matrix_complex_get_row(v: [*c]gsl_vector_complex, m: [*c]const gsl_matrix_complex, i: usize) c_int;
pub extern fn gsl_matrix_complex_get_col(v: [*c]gsl_vector_complex, m: [*c]const gsl_matrix_complex, j: usize) c_int;
pub extern fn gsl_matrix_complex_set_row(m: [*c]gsl_matrix_complex, i: usize, v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_matrix_complex_set_col(m: [*c]gsl_matrix_complex, j: usize, v: [*c]const gsl_vector_complex) c_int;
pub extern fn gsl_matrix_complex_get(m: [*c]const gsl_matrix_complex, i: usize, j: usize) gsl_complex;
pub extern fn gsl_matrix_complex_set(m: [*c]gsl_matrix_complex, i: usize, j: usize, x: gsl_complex) void;
pub extern fn gsl_matrix_complex_ptr(m: [*c]gsl_matrix_complex, i: usize, j: usize) [*c]gsl_complex;
pub extern fn gsl_matrix_complex_const_ptr(m: [*c]const gsl_matrix_complex, i: usize, j: usize) [*c]const gsl_complex;
pub const gsl_matrix_complex_float = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
    block: [*c]gsl_block_complex_float = @import("std").mem.zeroes([*c]gsl_block_complex_float),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_complex_float_view = extern struct {
    matrix: gsl_matrix_complex_float = @import("std").mem.zeroes(gsl_matrix_complex_float),
};
pub const gsl_matrix_complex_float_view = _gsl_matrix_complex_float_view;
pub const _gsl_matrix_complex_float_const_view = extern struct {
    matrix: gsl_matrix_complex_float = @import("std").mem.zeroes(gsl_matrix_complex_float),
};
pub const gsl_matrix_complex_float_const_view = _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_alloc(n1: usize, n2: usize) [*c]gsl_matrix_complex_float;
pub extern fn gsl_matrix_complex_float_calloc(n1: usize, n2: usize) [*c]gsl_matrix_complex_float;
pub extern fn gsl_matrix_complex_float_alloc_from_block(b: [*c]gsl_block_complex_float, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_complex_float;
pub extern fn gsl_matrix_complex_float_alloc_from_matrix(b: [*c]gsl_matrix_complex_float, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_complex_float;
pub extern fn gsl_vector_complex_float_alloc_row_from_matrix(m: [*c]gsl_matrix_complex_float, i: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_vector_complex_float_alloc_col_from_matrix(m: [*c]gsl_matrix_complex_float, j: usize) [*c]gsl_vector_complex_float;
pub extern fn gsl_matrix_complex_float_free(m: [*c]gsl_matrix_complex_float) void;
pub extern fn gsl_matrix_complex_float_submatrix(m: [*c]gsl_matrix_complex_float, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_float_view;
pub extern fn gsl_matrix_complex_float_row(m: [*c]gsl_matrix_complex_float, i: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_column(m: [*c]gsl_matrix_complex_float, j: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_diagonal(m: [*c]gsl_matrix_complex_float) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_subdiagonal(m: [*c]gsl_matrix_complex_float, k: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_superdiagonal(m: [*c]gsl_matrix_complex_float, k: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_subrow(m: [*c]gsl_matrix_complex_float, i: usize, offset: usize, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_subcolumn(m: [*c]gsl_matrix_complex_float, j: usize, offset: usize, n: usize) _gsl_vector_complex_float_view;
pub extern fn gsl_matrix_complex_float_view_array(base: [*c]f32, n1: usize, n2: usize) _gsl_matrix_complex_float_view;
pub extern fn gsl_matrix_complex_float_view_array_with_tda(base: [*c]f32, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_float_view;
pub extern fn gsl_matrix_complex_float_view_vector(v: [*c]gsl_vector_complex_float, n1: usize, n2: usize) _gsl_matrix_complex_float_view;
pub extern fn gsl_matrix_complex_float_view_vector_with_tda(v: [*c]gsl_vector_complex_float, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_float_view;
pub extern fn gsl_matrix_complex_float_const_submatrix(m: [*c]const gsl_matrix_complex_float, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_row(m: [*c]const gsl_matrix_complex_float, i: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_column(m: [*c]const gsl_matrix_complex_float, j: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_diagonal(m: [*c]const gsl_matrix_complex_float) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_subdiagonal(m: [*c]const gsl_matrix_complex_float, k: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_superdiagonal(m: [*c]const gsl_matrix_complex_float, k: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_subrow(m: [*c]const gsl_matrix_complex_float, i: usize, offset: usize, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_subcolumn(m: [*c]const gsl_matrix_complex_float, j: usize, offset: usize, n: usize) _gsl_vector_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_view_array(base: [*c]const f32, n1: usize, n2: usize) _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_view_array_with_tda(base: [*c]const f32, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_view_vector(v: [*c]const gsl_vector_complex_float, n1: usize, n2: usize) _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_const_view_vector_with_tda(v: [*c]const gsl_vector_complex_float, n1: usize, n2: usize, tda: usize) _gsl_matrix_complex_float_const_view;
pub extern fn gsl_matrix_complex_float_set_zero(m: [*c]gsl_matrix_complex_float) void;
pub extern fn gsl_matrix_complex_float_set_identity(m: [*c]gsl_matrix_complex_float) void;
pub extern fn gsl_matrix_complex_float_set_all(m: [*c]gsl_matrix_complex_float, x: gsl_complex_float) void;
pub extern fn gsl_matrix_complex_float_fread(stream: ?*FILE, m: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_complex_float, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_complex_float_memcpy(dest: [*c]gsl_matrix_complex_float, src: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_swap(m1: [*c]gsl_matrix_complex_float, m2: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex_float, src: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_swap_rows(m: [*c]gsl_matrix_complex_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_float_swap_columns(m: [*c]gsl_matrix_complex_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_float_swap_rowcol(m: [*c]gsl_matrix_complex_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_complex_float_transpose(m: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_transpose_memcpy(dest: [*c]gsl_matrix_complex_float, src: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_complex_float, src: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_conjtrans_memcpy(dest: [*c]gsl_matrix_complex_float, src: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_equal(a: [*c]const gsl_matrix_complex_float, b: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_isnull(m: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_ispos(m: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_isneg(m: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_isnonneg(m: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_add(a: [*c]gsl_matrix_complex_float, b: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_sub(a: [*c]gsl_matrix_complex_float, b: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_mul_elements(a: [*c]gsl_matrix_complex_float, b: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_div_elements(a: [*c]gsl_matrix_complex_float, b: [*c]const gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_scale(a: [*c]gsl_matrix_complex_float, x: gsl_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_scale_rows(a: [*c]gsl_matrix_complex_float, x: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_scale_columns(a: [*c]gsl_matrix_complex_float, x: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_add_constant(a: [*c]gsl_matrix_complex_float, x: gsl_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_add_diagonal(a: [*c]gsl_matrix_complex_float, x: gsl_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_conjugate(a: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_get_row(v: [*c]gsl_vector_complex_float, m: [*c]const gsl_matrix_complex_float, i: usize) c_int;
pub extern fn gsl_matrix_complex_float_get_col(v: [*c]gsl_vector_complex_float, m: [*c]const gsl_matrix_complex_float, j: usize) c_int;
pub extern fn gsl_matrix_complex_float_set_row(m: [*c]gsl_matrix_complex_float, i: usize, v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_set_col(m: [*c]gsl_matrix_complex_float, j: usize, v: [*c]const gsl_vector_complex_float) c_int;
pub extern fn gsl_matrix_complex_float_get(m: [*c]const gsl_matrix_complex_float, i: usize, j: usize) gsl_complex_float;
pub extern fn gsl_matrix_complex_float_set(m: [*c]gsl_matrix_complex_float, i: usize, j: usize, x: gsl_complex_float) void;
pub extern fn gsl_matrix_complex_float_ptr(m: [*c]gsl_matrix_complex_float, i: usize, j: usize) [*c]gsl_complex_float;
pub extern fn gsl_matrix_complex_float_const_ptr(m: [*c]const gsl_matrix_complex_float, i: usize, j: usize) [*c]const gsl_complex_float;
pub const gsl_matrix_long_double = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_longdouble = @import("std").mem.zeroes([*c]c_longdouble),
    block: [*c]gsl_block_long_double = @import("std").mem.zeroes([*c]gsl_block_long_double),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_long_double_view = extern struct {
    matrix: gsl_matrix_long_double = @import("std").mem.zeroes(gsl_matrix_long_double),
};
pub const gsl_matrix_long_double_view = _gsl_matrix_long_double_view;
pub const _gsl_matrix_long_double_const_view = extern struct {
    matrix: gsl_matrix_long_double = @import("std").mem.zeroes(gsl_matrix_long_double),
};
pub const gsl_matrix_long_double_const_view = _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_alloc(n1: usize, n2: usize) [*c]gsl_matrix_long_double;
pub extern fn gsl_matrix_long_double_calloc(n1: usize, n2: usize) [*c]gsl_matrix_long_double;
pub extern fn gsl_matrix_long_double_alloc_from_block(b: [*c]gsl_block_long_double, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_long_double;
pub extern fn gsl_matrix_long_double_alloc_from_matrix(m: [*c]gsl_matrix_long_double, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_long_double;
pub extern fn gsl_vector_long_double_alloc_row_from_matrix(m: [*c]gsl_matrix_long_double, i: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_vector_long_double_alloc_col_from_matrix(m: [*c]gsl_matrix_long_double, j: usize) [*c]gsl_vector_long_double;
pub extern fn gsl_matrix_long_double_free(m: [*c]gsl_matrix_long_double) void;
pub extern fn gsl_matrix_long_double_submatrix(m: [*c]gsl_matrix_long_double, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_long_double_view;
pub extern fn gsl_matrix_long_double_row(m: [*c]gsl_matrix_long_double, i: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_column(m: [*c]gsl_matrix_long_double, j: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_diagonal(m: [*c]gsl_matrix_long_double) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_subdiagonal(m: [*c]gsl_matrix_long_double, k: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_superdiagonal(m: [*c]gsl_matrix_long_double, k: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_subrow(m: [*c]gsl_matrix_long_double, i: usize, offset: usize, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_subcolumn(m: [*c]gsl_matrix_long_double, j: usize, offset: usize, n: usize) _gsl_vector_long_double_view;
pub extern fn gsl_matrix_long_double_view_array(base: [*c]c_longdouble, n1: usize, n2: usize) _gsl_matrix_long_double_view;
pub extern fn gsl_matrix_long_double_view_array_with_tda(base: [*c]c_longdouble, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_double_view;
pub extern fn gsl_matrix_long_double_view_vector(v: [*c]gsl_vector_long_double, n1: usize, n2: usize) _gsl_matrix_long_double_view;
pub extern fn gsl_matrix_long_double_view_vector_with_tda(v: [*c]gsl_vector_long_double, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_double_view;
pub extern fn gsl_matrix_long_double_const_submatrix(m: [*c]const gsl_matrix_long_double, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_row(m: [*c]const gsl_matrix_long_double, i: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_column(m: [*c]const gsl_matrix_long_double, j: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_diagonal(m: [*c]const gsl_matrix_long_double) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_subdiagonal(m: [*c]const gsl_matrix_long_double, k: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_superdiagonal(m: [*c]const gsl_matrix_long_double, k: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_subrow(m: [*c]const gsl_matrix_long_double, i: usize, offset: usize, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_subcolumn(m: [*c]const gsl_matrix_long_double, j: usize, offset: usize, n: usize) _gsl_vector_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_view_array(base: [*c]const c_longdouble, n1: usize, n2: usize) _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_view_array_with_tda(base: [*c]const c_longdouble, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_view_vector(v: [*c]const gsl_vector_long_double, n1: usize, n2: usize) _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_const_view_vector_with_tda(v: [*c]const gsl_vector_long_double, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_double_const_view;
pub extern fn gsl_matrix_long_double_set_zero(m: [*c]gsl_matrix_long_double) void;
pub extern fn gsl_matrix_long_double_set_identity(m: [*c]gsl_matrix_long_double) void;
pub extern fn gsl_matrix_long_double_set_all(m: [*c]gsl_matrix_long_double, x: c_longdouble) void;
pub extern fn gsl_matrix_long_double_fread(stream: ?*FILE, m: [*c]gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_long_double, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_long_double_memcpy(dest: [*c]gsl_matrix_long_double, src: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_swap(m1: [*c]gsl_matrix_long_double, m2: [*c]gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_long_double, src: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_swap_rows(m: [*c]gsl_matrix_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_double_swap_columns(m: [*c]gsl_matrix_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_double_swap_rowcol(m: [*c]gsl_matrix_long_double, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_double_transpose(m: [*c]gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_transpose_memcpy(dest: [*c]gsl_matrix_long_double, src: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_long_double, src: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_max(m: [*c]const gsl_matrix_long_double) c_longdouble;
pub extern fn gsl_matrix_long_double_min(m: [*c]const gsl_matrix_long_double) c_longdouble;
pub extern fn gsl_matrix_long_double_minmax(m: [*c]const gsl_matrix_long_double, min_out: [*c]c_longdouble, max_out: [*c]c_longdouble) void;
pub extern fn gsl_matrix_long_double_max_index(m: [*c]const gsl_matrix_long_double, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_long_double_min_index(m: [*c]const gsl_matrix_long_double, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_long_double_minmax_index(m: [*c]const gsl_matrix_long_double, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_long_double_equal(a: [*c]const gsl_matrix_long_double, b: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_isnull(m: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_ispos(m: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_isneg(m: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_isnonneg(m: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_norm1(m: [*c]const gsl_matrix_long_double) c_longdouble;
pub extern fn gsl_matrix_long_double_add(a: [*c]gsl_matrix_long_double, b: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_sub(a: [*c]gsl_matrix_long_double, b: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_mul_elements(a: [*c]gsl_matrix_long_double, b: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_div_elements(a: [*c]gsl_matrix_long_double, b: [*c]const gsl_matrix_long_double) c_int;
pub extern fn gsl_matrix_long_double_scale(a: [*c]gsl_matrix_long_double, x: c_longdouble) c_int;
pub extern fn gsl_matrix_long_double_scale_rows(a: [*c]gsl_matrix_long_double, x: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_matrix_long_double_scale_columns(a: [*c]gsl_matrix_long_double, x: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_matrix_long_double_add_constant(a: [*c]gsl_matrix_long_double, x: c_longdouble) c_int;
pub extern fn gsl_matrix_long_double_add_diagonal(a: [*c]gsl_matrix_long_double, x: c_longdouble) c_int;
pub extern fn gsl_matrix_long_double_get_row(v: [*c]gsl_vector_long_double, m: [*c]const gsl_matrix_long_double, i: usize) c_int;
pub extern fn gsl_matrix_long_double_get_col(v: [*c]gsl_vector_long_double, m: [*c]const gsl_matrix_long_double, j: usize) c_int;
pub extern fn gsl_matrix_long_double_set_row(m: [*c]gsl_matrix_long_double, i: usize, v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_matrix_long_double_set_col(m: [*c]gsl_matrix_long_double, j: usize, v: [*c]const gsl_vector_long_double) c_int;
pub extern fn gsl_matrix_long_double_get(m: [*c]const gsl_matrix_long_double, i: usize, j: usize) c_longdouble;
pub extern fn gsl_matrix_long_double_set(m: [*c]gsl_matrix_long_double, i: usize, j: usize, x: c_longdouble) void;
pub extern fn gsl_matrix_long_double_ptr(m: [*c]gsl_matrix_long_double, i: usize, j: usize) [*c]c_longdouble;
pub extern fn gsl_matrix_long_double_const_ptr(m: [*c]const gsl_matrix_long_double, i: usize, j: usize) [*c]const c_longdouble;
pub const gsl_matrix = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    block: [*c]gsl_block = @import("std").mem.zeroes([*c]gsl_block),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_view = extern struct {
    matrix: gsl_matrix = @import("std").mem.zeroes(gsl_matrix),
};
pub const gsl_matrix_view = _gsl_matrix_view;
pub const _gsl_matrix_const_view = extern struct {
    matrix: gsl_matrix = @import("std").mem.zeroes(gsl_matrix),
};
pub const gsl_matrix_const_view = _gsl_matrix_const_view;
pub extern fn gsl_matrix_alloc(n1: usize, n2: usize) [*c]gsl_matrix;
pub extern fn gsl_matrix_calloc(n1: usize, n2: usize) [*c]gsl_matrix;
pub extern fn gsl_matrix_alloc_from_block(b: [*c]gsl_block, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix;
pub extern fn gsl_matrix_alloc_from_matrix(m: [*c]gsl_matrix, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix;
pub extern fn gsl_vector_alloc_row_from_matrix(m: [*c]gsl_matrix, i: usize) [*c]gsl_vector;
pub extern fn gsl_vector_alloc_col_from_matrix(m: [*c]gsl_matrix, j: usize) [*c]gsl_vector;
pub extern fn gsl_matrix_free(m: [*c]gsl_matrix) void;
pub extern fn gsl_matrix_submatrix(m: [*c]gsl_matrix, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_view;
pub extern fn gsl_matrix_row(m: [*c]gsl_matrix, i: usize) _gsl_vector_view;
pub extern fn gsl_matrix_column(m: [*c]gsl_matrix, j: usize) _gsl_vector_view;
pub extern fn gsl_matrix_diagonal(m: [*c]gsl_matrix) _gsl_vector_view;
pub extern fn gsl_matrix_subdiagonal(m: [*c]gsl_matrix, k: usize) _gsl_vector_view;
pub extern fn gsl_matrix_superdiagonal(m: [*c]gsl_matrix, k: usize) _gsl_vector_view;
pub extern fn gsl_matrix_subrow(m: [*c]gsl_matrix, i: usize, offset: usize, n: usize) _gsl_vector_view;
pub extern fn gsl_matrix_subcolumn(m: [*c]gsl_matrix, j: usize, offset: usize, n: usize) _gsl_vector_view;
pub extern fn gsl_matrix_view_array(base: [*c]f64, n1: usize, n2: usize) _gsl_matrix_view;
pub extern fn gsl_matrix_view_array_with_tda(base: [*c]f64, n1: usize, n2: usize, tda: usize) _gsl_matrix_view;
pub extern fn gsl_matrix_view_vector(v: [*c]gsl_vector, n1: usize, n2: usize) _gsl_matrix_view;
pub extern fn gsl_matrix_view_vector_with_tda(v: [*c]gsl_vector, n1: usize, n2: usize, tda: usize) _gsl_matrix_view;
pub extern fn gsl_matrix_const_submatrix(m: [*c]const gsl_matrix, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_const_view;
pub extern fn gsl_matrix_const_row(m: [*c]const gsl_matrix, i: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_column(m: [*c]const gsl_matrix, j: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_diagonal(m: [*c]const gsl_matrix) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_subdiagonal(m: [*c]const gsl_matrix, k: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_superdiagonal(m: [*c]const gsl_matrix, k: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_subrow(m: [*c]const gsl_matrix, i: usize, offset: usize, n: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_subcolumn(m: [*c]const gsl_matrix, j: usize, offset: usize, n: usize) _gsl_vector_const_view;
pub extern fn gsl_matrix_const_view_array(base: [*c]const f64, n1: usize, n2: usize) _gsl_matrix_const_view;
pub extern fn gsl_matrix_const_view_array_with_tda(base: [*c]const f64, n1: usize, n2: usize, tda: usize) _gsl_matrix_const_view;
pub extern fn gsl_matrix_const_view_vector(v: [*c]const gsl_vector, n1: usize, n2: usize) _gsl_matrix_const_view;
pub extern fn gsl_matrix_const_view_vector_with_tda(v: [*c]const gsl_vector, n1: usize, n2: usize, tda: usize) _gsl_matrix_const_view;
pub extern fn gsl_matrix_set_zero(m: [*c]gsl_matrix) void;
pub extern fn gsl_matrix_set_identity(m: [*c]gsl_matrix) void;
pub extern fn gsl_matrix_set_all(m: [*c]gsl_matrix, x: f64) void;
pub extern fn gsl_matrix_fread(stream: ?*FILE, m: [*c]gsl_matrix) c_int;
pub extern fn gsl_matrix_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_fscanf(stream: ?*FILE, m: [*c]gsl_matrix) c_int;
pub extern fn gsl_matrix_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_memcpy(dest: [*c]gsl_matrix, src: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_swap(m1: [*c]gsl_matrix, m2: [*c]gsl_matrix) c_int;
pub extern fn gsl_matrix_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix, src: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_swap_rows(m: [*c]gsl_matrix, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_swap_columns(m: [*c]gsl_matrix, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_swap_rowcol(m: [*c]gsl_matrix, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_transpose(m: [*c]gsl_matrix) c_int;
pub extern fn gsl_matrix_transpose_memcpy(dest: [*c]gsl_matrix, src: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix, src: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_max(m: [*c]const gsl_matrix) f64;
pub extern fn gsl_matrix_min(m: [*c]const gsl_matrix) f64;
pub extern fn gsl_matrix_minmax(m: [*c]const gsl_matrix, min_out: [*c]f64, max_out: [*c]f64) void;
pub extern fn gsl_matrix_max_index(m: [*c]const gsl_matrix, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_min_index(m: [*c]const gsl_matrix, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_minmax_index(m: [*c]const gsl_matrix, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_equal(a: [*c]const gsl_matrix, b: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_isnull(m: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_ispos(m: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_isneg(m: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_isnonneg(m: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_norm1(m: [*c]const gsl_matrix) f64;
pub extern fn gsl_matrix_add(a: [*c]gsl_matrix, b: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_sub(a: [*c]gsl_matrix, b: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_mul_elements(a: [*c]gsl_matrix, b: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_div_elements(a: [*c]gsl_matrix, b: [*c]const gsl_matrix) c_int;
pub extern fn gsl_matrix_scale(a: [*c]gsl_matrix, x: f64) c_int;
pub extern fn gsl_matrix_scale_rows(a: [*c]gsl_matrix, x: [*c]const gsl_vector) c_int;
pub extern fn gsl_matrix_scale_columns(a: [*c]gsl_matrix, x: [*c]const gsl_vector) c_int;
pub extern fn gsl_matrix_add_constant(a: [*c]gsl_matrix, x: f64) c_int;
pub extern fn gsl_matrix_add_diagonal(a: [*c]gsl_matrix, x: f64) c_int;
pub extern fn gsl_matrix_get_row(v: [*c]gsl_vector, m: [*c]const gsl_matrix, i: usize) c_int;
pub extern fn gsl_matrix_get_col(v: [*c]gsl_vector, m: [*c]const gsl_matrix, j: usize) c_int;
pub extern fn gsl_matrix_set_row(m: [*c]gsl_matrix, i: usize, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_matrix_set_col(m: [*c]gsl_matrix, j: usize, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_matrix_get(m: [*c]const gsl_matrix, i: usize, j: usize) f64;
pub extern fn gsl_matrix_set(m: [*c]gsl_matrix, i: usize, j: usize, x: f64) void;
pub extern fn gsl_matrix_ptr(m: [*c]gsl_matrix, i: usize, j: usize) [*c]f64;
pub extern fn gsl_matrix_const_ptr(m: [*c]const gsl_matrix, i: usize, j: usize) [*c]const f64;
pub const gsl_matrix_float = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]f32 = @import("std").mem.zeroes([*c]f32),
    block: [*c]gsl_block_float = @import("std").mem.zeroes([*c]gsl_block_float),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_float_view = extern struct {
    matrix: gsl_matrix_float = @import("std").mem.zeroes(gsl_matrix_float),
};
pub const gsl_matrix_float_view = _gsl_matrix_float_view;
pub const _gsl_matrix_float_const_view = extern struct {
    matrix: gsl_matrix_float = @import("std").mem.zeroes(gsl_matrix_float),
};
pub const gsl_matrix_float_const_view = _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_alloc(n1: usize, n2: usize) [*c]gsl_matrix_float;
pub extern fn gsl_matrix_float_calloc(n1: usize, n2: usize) [*c]gsl_matrix_float;
pub extern fn gsl_matrix_float_alloc_from_block(b: [*c]gsl_block_float, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_float;
pub extern fn gsl_matrix_float_alloc_from_matrix(m: [*c]gsl_matrix_float, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_float;
pub extern fn gsl_vector_float_alloc_row_from_matrix(m: [*c]gsl_matrix_float, i: usize) [*c]gsl_vector_float;
pub extern fn gsl_vector_float_alloc_col_from_matrix(m: [*c]gsl_matrix_float, j: usize) [*c]gsl_vector_float;
pub extern fn gsl_matrix_float_free(m: [*c]gsl_matrix_float) void;
pub extern fn gsl_matrix_float_submatrix(m: [*c]gsl_matrix_float, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_float_view;
pub extern fn gsl_matrix_float_row(m: [*c]gsl_matrix_float, i: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_column(m: [*c]gsl_matrix_float, j: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_diagonal(m: [*c]gsl_matrix_float) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_subdiagonal(m: [*c]gsl_matrix_float, k: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_superdiagonal(m: [*c]gsl_matrix_float, k: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_subrow(m: [*c]gsl_matrix_float, i: usize, offset: usize, n: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_subcolumn(m: [*c]gsl_matrix_float, j: usize, offset: usize, n: usize) _gsl_vector_float_view;
pub extern fn gsl_matrix_float_view_array(base: [*c]f32, n1: usize, n2: usize) _gsl_matrix_float_view;
pub extern fn gsl_matrix_float_view_array_with_tda(base: [*c]f32, n1: usize, n2: usize, tda: usize) _gsl_matrix_float_view;
pub extern fn gsl_matrix_float_view_vector(v: [*c]gsl_vector_float, n1: usize, n2: usize) _gsl_matrix_float_view;
pub extern fn gsl_matrix_float_view_vector_with_tda(v: [*c]gsl_vector_float, n1: usize, n2: usize, tda: usize) _gsl_matrix_float_view;
pub extern fn gsl_matrix_float_const_submatrix(m: [*c]const gsl_matrix_float, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_const_row(m: [*c]const gsl_matrix_float, i: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_column(m: [*c]const gsl_matrix_float, j: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_diagonal(m: [*c]const gsl_matrix_float) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_subdiagonal(m: [*c]const gsl_matrix_float, k: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_superdiagonal(m: [*c]const gsl_matrix_float, k: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_subrow(m: [*c]const gsl_matrix_float, i: usize, offset: usize, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_subcolumn(m: [*c]const gsl_matrix_float, j: usize, offset: usize, n: usize) _gsl_vector_float_const_view;
pub extern fn gsl_matrix_float_const_view_array(base: [*c]const f32, n1: usize, n2: usize) _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_const_view_array_with_tda(base: [*c]const f32, n1: usize, n2: usize, tda: usize) _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_const_view_vector(v: [*c]const gsl_vector_float, n1: usize, n2: usize) _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_const_view_vector_with_tda(v: [*c]const gsl_vector_float, n1: usize, n2: usize, tda: usize) _gsl_matrix_float_const_view;
pub extern fn gsl_matrix_float_set_zero(m: [*c]gsl_matrix_float) void;
pub extern fn gsl_matrix_float_set_identity(m: [*c]gsl_matrix_float) void;
pub extern fn gsl_matrix_float_set_all(m: [*c]gsl_matrix_float, x: f32) void;
pub extern fn gsl_matrix_float_fread(stream: ?*FILE, m: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_float, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_float_memcpy(dest: [*c]gsl_matrix_float, src: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_swap(m1: [*c]gsl_matrix_float, m2: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_float, src: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_swap_rows(m: [*c]gsl_matrix_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_float_swap_columns(m: [*c]gsl_matrix_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_float_swap_rowcol(m: [*c]gsl_matrix_float, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_float_transpose(m: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_transpose_memcpy(dest: [*c]gsl_matrix_float, src: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_float, src: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_max(m: [*c]const gsl_matrix_float) f32;
pub extern fn gsl_matrix_float_min(m: [*c]const gsl_matrix_float) f32;
pub extern fn gsl_matrix_float_minmax(m: [*c]const gsl_matrix_float, min_out: [*c]f32, max_out: [*c]f32) void;
pub extern fn gsl_matrix_float_max_index(m: [*c]const gsl_matrix_float, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_float_min_index(m: [*c]const gsl_matrix_float, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_float_minmax_index(m: [*c]const gsl_matrix_float, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_float_equal(a: [*c]const gsl_matrix_float, b: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_isnull(m: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_ispos(m: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_isneg(m: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_isnonneg(m: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_norm1(m: [*c]const gsl_matrix_float) f32;
pub extern fn gsl_matrix_float_add(a: [*c]gsl_matrix_float, b: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_sub(a: [*c]gsl_matrix_float, b: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_mul_elements(a: [*c]gsl_matrix_float, b: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_div_elements(a: [*c]gsl_matrix_float, b: [*c]const gsl_matrix_float) c_int;
pub extern fn gsl_matrix_float_scale(a: [*c]gsl_matrix_float, x: f32) c_int;
pub extern fn gsl_matrix_float_scale_rows(a: [*c]gsl_matrix_float, x: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_matrix_float_scale_columns(a: [*c]gsl_matrix_float, x: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_matrix_float_add_constant(a: [*c]gsl_matrix_float, x: f32) c_int;
pub extern fn gsl_matrix_float_add_diagonal(a: [*c]gsl_matrix_float, x: f32) c_int;
pub extern fn gsl_matrix_float_get_row(v: [*c]gsl_vector_float, m: [*c]const gsl_matrix_float, i: usize) c_int;
pub extern fn gsl_matrix_float_get_col(v: [*c]gsl_vector_float, m: [*c]const gsl_matrix_float, j: usize) c_int;
pub extern fn gsl_matrix_float_set_row(m: [*c]gsl_matrix_float, i: usize, v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_matrix_float_set_col(m: [*c]gsl_matrix_float, j: usize, v: [*c]const gsl_vector_float) c_int;
pub extern fn gsl_matrix_float_get(m: [*c]const gsl_matrix_float, i: usize, j: usize) f32;
pub extern fn gsl_matrix_float_set(m: [*c]gsl_matrix_float, i: usize, j: usize, x: f32) void;
pub extern fn gsl_matrix_float_ptr(m: [*c]gsl_matrix_float, i: usize, j: usize) [*c]f32;
pub extern fn gsl_matrix_float_const_ptr(m: [*c]const gsl_matrix_float, i: usize, j: usize) [*c]const f32;
pub const gsl_matrix_ulong = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ulong = @import("std").mem.zeroes([*c]c_ulong),
    block: [*c]gsl_block_ulong = @import("std").mem.zeroes([*c]gsl_block_ulong),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_ulong_view = extern struct {
    matrix: gsl_matrix_ulong = @import("std").mem.zeroes(gsl_matrix_ulong),
};
pub const gsl_matrix_ulong_view = _gsl_matrix_ulong_view;
pub const _gsl_matrix_ulong_const_view = extern struct {
    matrix: gsl_matrix_ulong = @import("std").mem.zeroes(gsl_matrix_ulong),
};
pub const gsl_matrix_ulong_const_view = _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_alloc(n1: usize, n2: usize) [*c]gsl_matrix_ulong;
pub extern fn gsl_matrix_ulong_calloc(n1: usize, n2: usize) [*c]gsl_matrix_ulong;
pub extern fn gsl_matrix_ulong_alloc_from_block(b: [*c]gsl_block_ulong, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_ulong;
pub extern fn gsl_matrix_ulong_alloc_from_matrix(m: [*c]gsl_matrix_ulong, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_ulong;
pub extern fn gsl_vector_ulong_alloc_row_from_matrix(m: [*c]gsl_matrix_ulong, i: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_vector_ulong_alloc_col_from_matrix(m: [*c]gsl_matrix_ulong, j: usize) [*c]gsl_vector_ulong;
pub extern fn gsl_matrix_ulong_free(m: [*c]gsl_matrix_ulong) void;
pub extern fn gsl_matrix_ulong_submatrix(m: [*c]gsl_matrix_ulong, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_ulong_view;
pub extern fn gsl_matrix_ulong_row(m: [*c]gsl_matrix_ulong, i: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_column(m: [*c]gsl_matrix_ulong, j: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_diagonal(m: [*c]gsl_matrix_ulong) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_subdiagonal(m: [*c]gsl_matrix_ulong, k: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_superdiagonal(m: [*c]gsl_matrix_ulong, k: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_subrow(m: [*c]gsl_matrix_ulong, i: usize, offset: usize, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_subcolumn(m: [*c]gsl_matrix_ulong, j: usize, offset: usize, n: usize) _gsl_vector_ulong_view;
pub extern fn gsl_matrix_ulong_view_array(base: [*c]c_ulong, n1: usize, n2: usize) _gsl_matrix_ulong_view;
pub extern fn gsl_matrix_ulong_view_array_with_tda(base: [*c]c_ulong, n1: usize, n2: usize, tda: usize) _gsl_matrix_ulong_view;
pub extern fn gsl_matrix_ulong_view_vector(v: [*c]gsl_vector_ulong, n1: usize, n2: usize) _gsl_matrix_ulong_view;
pub extern fn gsl_matrix_ulong_view_vector_with_tda(v: [*c]gsl_vector_ulong, n1: usize, n2: usize, tda: usize) _gsl_matrix_ulong_view;
pub extern fn gsl_matrix_ulong_const_submatrix(m: [*c]const gsl_matrix_ulong, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_row(m: [*c]const gsl_matrix_ulong, i: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_column(m: [*c]const gsl_matrix_ulong, j: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_diagonal(m: [*c]const gsl_matrix_ulong) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_subdiagonal(m: [*c]const gsl_matrix_ulong, k: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_superdiagonal(m: [*c]const gsl_matrix_ulong, k: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_subrow(m: [*c]const gsl_matrix_ulong, i: usize, offset: usize, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_subcolumn(m: [*c]const gsl_matrix_ulong, j: usize, offset: usize, n: usize) _gsl_vector_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_view_array(base: [*c]const c_ulong, n1: usize, n2: usize) _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_view_array_with_tda(base: [*c]const c_ulong, n1: usize, n2: usize, tda: usize) _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_view_vector(v: [*c]const gsl_vector_ulong, n1: usize, n2: usize) _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_const_view_vector_with_tda(v: [*c]const gsl_vector_ulong, n1: usize, n2: usize, tda: usize) _gsl_matrix_ulong_const_view;
pub extern fn gsl_matrix_ulong_set_zero(m: [*c]gsl_matrix_ulong) void;
pub extern fn gsl_matrix_ulong_set_identity(m: [*c]gsl_matrix_ulong) void;
pub extern fn gsl_matrix_ulong_set_all(m: [*c]gsl_matrix_ulong, x: c_ulong) void;
pub extern fn gsl_matrix_ulong_fread(stream: ?*FILE, m: [*c]gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_ulong, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_ulong_memcpy(dest: [*c]gsl_matrix_ulong, src: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_swap(m1: [*c]gsl_matrix_ulong, m2: [*c]gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_ulong, src: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_swap_rows(m: [*c]gsl_matrix_ulong, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ulong_swap_columns(m: [*c]gsl_matrix_ulong, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ulong_swap_rowcol(m: [*c]gsl_matrix_ulong, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ulong_transpose(m: [*c]gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_transpose_memcpy(dest: [*c]gsl_matrix_ulong, src: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_ulong, src: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_max(m: [*c]const gsl_matrix_ulong) c_ulong;
pub extern fn gsl_matrix_ulong_min(m: [*c]const gsl_matrix_ulong) c_ulong;
pub extern fn gsl_matrix_ulong_minmax(m: [*c]const gsl_matrix_ulong, min_out: [*c]c_ulong, max_out: [*c]c_ulong) void;
pub extern fn gsl_matrix_ulong_max_index(m: [*c]const gsl_matrix_ulong, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_ulong_min_index(m: [*c]const gsl_matrix_ulong, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_ulong_minmax_index(m: [*c]const gsl_matrix_ulong, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_ulong_equal(a: [*c]const gsl_matrix_ulong, b: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_isnull(m: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_ispos(m: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_isneg(m: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_isnonneg(m: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_norm1(m: [*c]const gsl_matrix_ulong) c_ulong;
pub extern fn gsl_matrix_ulong_add(a: [*c]gsl_matrix_ulong, b: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_sub(a: [*c]gsl_matrix_ulong, b: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_mul_elements(a: [*c]gsl_matrix_ulong, b: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_div_elements(a: [*c]gsl_matrix_ulong, b: [*c]const gsl_matrix_ulong) c_int;
pub extern fn gsl_matrix_ulong_scale(a: [*c]gsl_matrix_ulong, x: c_ulong) c_int;
pub extern fn gsl_matrix_ulong_scale_rows(a: [*c]gsl_matrix_ulong, x: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_matrix_ulong_scale_columns(a: [*c]gsl_matrix_ulong, x: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_matrix_ulong_add_constant(a: [*c]gsl_matrix_ulong, x: c_ulong) c_int;
pub extern fn gsl_matrix_ulong_add_diagonal(a: [*c]gsl_matrix_ulong, x: c_ulong) c_int;
pub extern fn gsl_matrix_ulong_get_row(v: [*c]gsl_vector_ulong, m: [*c]const gsl_matrix_ulong, i: usize) c_int;
pub extern fn gsl_matrix_ulong_get_col(v: [*c]gsl_vector_ulong, m: [*c]const gsl_matrix_ulong, j: usize) c_int;
pub extern fn gsl_matrix_ulong_set_row(m: [*c]gsl_matrix_ulong, i: usize, v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_matrix_ulong_set_col(m: [*c]gsl_matrix_ulong, j: usize, v: [*c]const gsl_vector_ulong) c_int;
pub extern fn gsl_matrix_ulong_get(m: [*c]const gsl_matrix_ulong, i: usize, j: usize) c_ulong;
pub extern fn gsl_matrix_ulong_set(m: [*c]gsl_matrix_ulong, i: usize, j: usize, x: c_ulong) void;
pub extern fn gsl_matrix_ulong_ptr(m: [*c]gsl_matrix_ulong, i: usize, j: usize) [*c]c_ulong;
pub extern fn gsl_matrix_ulong_const_ptr(m: [*c]const gsl_matrix_ulong, i: usize, j: usize) [*c]const c_ulong;
pub const gsl_matrix_long = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_long = @import("std").mem.zeroes([*c]c_long),
    block: [*c]gsl_block_long = @import("std").mem.zeroes([*c]gsl_block_long),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_long_view = extern struct {
    matrix: gsl_matrix_long = @import("std").mem.zeroes(gsl_matrix_long),
};
pub const gsl_matrix_long_view = _gsl_matrix_long_view;
pub const _gsl_matrix_long_const_view = extern struct {
    matrix: gsl_matrix_long = @import("std").mem.zeroes(gsl_matrix_long),
};
pub const gsl_matrix_long_const_view = _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_alloc(n1: usize, n2: usize) [*c]gsl_matrix_long;
pub extern fn gsl_matrix_long_calloc(n1: usize, n2: usize) [*c]gsl_matrix_long;
pub extern fn gsl_matrix_long_alloc_from_block(b: [*c]gsl_block_long, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_long;
pub extern fn gsl_matrix_long_alloc_from_matrix(m: [*c]gsl_matrix_long, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_long;
pub extern fn gsl_vector_long_alloc_row_from_matrix(m: [*c]gsl_matrix_long, i: usize) [*c]gsl_vector_long;
pub extern fn gsl_vector_long_alloc_col_from_matrix(m: [*c]gsl_matrix_long, j: usize) [*c]gsl_vector_long;
pub extern fn gsl_matrix_long_free(m: [*c]gsl_matrix_long) void;
pub extern fn gsl_matrix_long_submatrix(m: [*c]gsl_matrix_long, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_long_view;
pub extern fn gsl_matrix_long_row(m: [*c]gsl_matrix_long, i: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_column(m: [*c]gsl_matrix_long, j: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_diagonal(m: [*c]gsl_matrix_long) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_subdiagonal(m: [*c]gsl_matrix_long, k: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_superdiagonal(m: [*c]gsl_matrix_long, k: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_subrow(m: [*c]gsl_matrix_long, i: usize, offset: usize, n: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_subcolumn(m: [*c]gsl_matrix_long, j: usize, offset: usize, n: usize) _gsl_vector_long_view;
pub extern fn gsl_matrix_long_view_array(base: [*c]c_long, n1: usize, n2: usize) _gsl_matrix_long_view;
pub extern fn gsl_matrix_long_view_array_with_tda(base: [*c]c_long, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_view;
pub extern fn gsl_matrix_long_view_vector(v: [*c]gsl_vector_long, n1: usize, n2: usize) _gsl_matrix_long_view;
pub extern fn gsl_matrix_long_view_vector_with_tda(v: [*c]gsl_vector_long, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_view;
pub extern fn gsl_matrix_long_const_submatrix(m: [*c]const gsl_matrix_long, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_const_row(m: [*c]const gsl_matrix_long, i: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_column(m: [*c]const gsl_matrix_long, j: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_diagonal(m: [*c]const gsl_matrix_long) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_subdiagonal(m: [*c]const gsl_matrix_long, k: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_superdiagonal(m: [*c]const gsl_matrix_long, k: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_subrow(m: [*c]const gsl_matrix_long, i: usize, offset: usize, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_subcolumn(m: [*c]const gsl_matrix_long, j: usize, offset: usize, n: usize) _gsl_vector_long_const_view;
pub extern fn gsl_matrix_long_const_view_array(base: [*c]const c_long, n1: usize, n2: usize) _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_const_view_array_with_tda(base: [*c]const c_long, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_const_view_vector(v: [*c]const gsl_vector_long, n1: usize, n2: usize) _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_const_view_vector_with_tda(v: [*c]const gsl_vector_long, n1: usize, n2: usize, tda: usize) _gsl_matrix_long_const_view;
pub extern fn gsl_matrix_long_set_zero(m: [*c]gsl_matrix_long) void;
pub extern fn gsl_matrix_long_set_identity(m: [*c]gsl_matrix_long) void;
pub extern fn gsl_matrix_long_set_all(m: [*c]gsl_matrix_long, x: c_long) void;
pub extern fn gsl_matrix_long_fread(stream: ?*FILE, m: [*c]gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_long, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_long_memcpy(dest: [*c]gsl_matrix_long, src: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_swap(m1: [*c]gsl_matrix_long, m2: [*c]gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_long, src: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_swap_rows(m: [*c]gsl_matrix_long, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_swap_columns(m: [*c]gsl_matrix_long, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_swap_rowcol(m: [*c]gsl_matrix_long, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_long_transpose(m: [*c]gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_transpose_memcpy(dest: [*c]gsl_matrix_long, src: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_long, src: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_max(m: [*c]const gsl_matrix_long) c_long;
pub extern fn gsl_matrix_long_min(m: [*c]const gsl_matrix_long) c_long;
pub extern fn gsl_matrix_long_minmax(m: [*c]const gsl_matrix_long, min_out: [*c]c_long, max_out: [*c]c_long) void;
pub extern fn gsl_matrix_long_max_index(m: [*c]const gsl_matrix_long, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_long_min_index(m: [*c]const gsl_matrix_long, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_long_minmax_index(m: [*c]const gsl_matrix_long, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_long_equal(a: [*c]const gsl_matrix_long, b: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_isnull(m: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_ispos(m: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_isneg(m: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_isnonneg(m: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_norm1(m: [*c]const gsl_matrix_long) c_long;
pub extern fn gsl_matrix_long_add(a: [*c]gsl_matrix_long, b: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_sub(a: [*c]gsl_matrix_long, b: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_mul_elements(a: [*c]gsl_matrix_long, b: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_div_elements(a: [*c]gsl_matrix_long, b: [*c]const gsl_matrix_long) c_int;
pub extern fn gsl_matrix_long_scale(a: [*c]gsl_matrix_long, x: c_long) c_int;
pub extern fn gsl_matrix_long_scale_rows(a: [*c]gsl_matrix_long, x: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_matrix_long_scale_columns(a: [*c]gsl_matrix_long, x: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_matrix_long_add_constant(a: [*c]gsl_matrix_long, x: c_long) c_int;
pub extern fn gsl_matrix_long_add_diagonal(a: [*c]gsl_matrix_long, x: c_long) c_int;
pub extern fn gsl_matrix_long_get_row(v: [*c]gsl_vector_long, m: [*c]const gsl_matrix_long, i: usize) c_int;
pub extern fn gsl_matrix_long_get_col(v: [*c]gsl_vector_long, m: [*c]const gsl_matrix_long, j: usize) c_int;
pub extern fn gsl_matrix_long_set_row(m: [*c]gsl_matrix_long, i: usize, v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_matrix_long_set_col(m: [*c]gsl_matrix_long, j: usize, v: [*c]const gsl_vector_long) c_int;
pub extern fn gsl_matrix_long_get(m: [*c]const gsl_matrix_long, i: usize, j: usize) c_long;
pub extern fn gsl_matrix_long_set(m: [*c]gsl_matrix_long, i: usize, j: usize, x: c_long) void;
pub extern fn gsl_matrix_long_ptr(m: [*c]gsl_matrix_long, i: usize, j: usize) [*c]c_long;
pub extern fn gsl_matrix_long_const_ptr(m: [*c]const gsl_matrix_long, i: usize, j: usize) [*c]const c_long;
pub const gsl_matrix_uint = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_uint = @import("std").mem.zeroes([*c]c_uint),
    block: [*c]gsl_block_uint = @import("std").mem.zeroes([*c]gsl_block_uint),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_uint_view = extern struct {
    matrix: gsl_matrix_uint = @import("std").mem.zeroes(gsl_matrix_uint),
};
pub const gsl_matrix_uint_view = _gsl_matrix_uint_view;
pub const _gsl_matrix_uint_const_view = extern struct {
    matrix: gsl_matrix_uint = @import("std").mem.zeroes(gsl_matrix_uint),
};
pub const gsl_matrix_uint_const_view = _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_alloc(n1: usize, n2: usize) [*c]gsl_matrix_uint;
pub extern fn gsl_matrix_uint_calloc(n1: usize, n2: usize) [*c]gsl_matrix_uint;
pub extern fn gsl_matrix_uint_alloc_from_block(b: [*c]gsl_block_uint, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_uint;
pub extern fn gsl_matrix_uint_alloc_from_matrix(m: [*c]gsl_matrix_uint, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_uint;
pub extern fn gsl_vector_uint_alloc_row_from_matrix(m: [*c]gsl_matrix_uint, i: usize) [*c]gsl_vector_uint;
pub extern fn gsl_vector_uint_alloc_col_from_matrix(m: [*c]gsl_matrix_uint, j: usize) [*c]gsl_vector_uint;
pub extern fn gsl_matrix_uint_free(m: [*c]gsl_matrix_uint) void;
pub extern fn gsl_matrix_uint_submatrix(m: [*c]gsl_matrix_uint, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_uint_view;
pub extern fn gsl_matrix_uint_row(m: [*c]gsl_matrix_uint, i: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_column(m: [*c]gsl_matrix_uint, j: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_diagonal(m: [*c]gsl_matrix_uint) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_subdiagonal(m: [*c]gsl_matrix_uint, k: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_superdiagonal(m: [*c]gsl_matrix_uint, k: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_subrow(m: [*c]gsl_matrix_uint, i: usize, offset: usize, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_subcolumn(m: [*c]gsl_matrix_uint, j: usize, offset: usize, n: usize) _gsl_vector_uint_view;
pub extern fn gsl_matrix_uint_view_array(base: [*c]c_uint, n1: usize, n2: usize) _gsl_matrix_uint_view;
pub extern fn gsl_matrix_uint_view_array_with_tda(base: [*c]c_uint, n1: usize, n2: usize, tda: usize) _gsl_matrix_uint_view;
pub extern fn gsl_matrix_uint_view_vector(v: [*c]gsl_vector_uint, n1: usize, n2: usize) _gsl_matrix_uint_view;
pub extern fn gsl_matrix_uint_view_vector_with_tda(v: [*c]gsl_vector_uint, n1: usize, n2: usize, tda: usize) _gsl_matrix_uint_view;
pub extern fn gsl_matrix_uint_const_submatrix(m: [*c]const gsl_matrix_uint, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_const_row(m: [*c]const gsl_matrix_uint, i: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_column(m: [*c]const gsl_matrix_uint, j: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_diagonal(m: [*c]const gsl_matrix_uint) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_subdiagonal(m: [*c]const gsl_matrix_uint, k: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_superdiagonal(m: [*c]const gsl_matrix_uint, k: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_subrow(m: [*c]const gsl_matrix_uint, i: usize, offset: usize, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_subcolumn(m: [*c]const gsl_matrix_uint, j: usize, offset: usize, n: usize) _gsl_vector_uint_const_view;
pub extern fn gsl_matrix_uint_const_view_array(base: [*c]const c_uint, n1: usize, n2: usize) _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_const_view_array_with_tda(base: [*c]const c_uint, n1: usize, n2: usize, tda: usize) _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_const_view_vector(v: [*c]const gsl_vector_uint, n1: usize, n2: usize) _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_const_view_vector_with_tda(v: [*c]const gsl_vector_uint, n1: usize, n2: usize, tda: usize) _gsl_matrix_uint_const_view;
pub extern fn gsl_matrix_uint_set_zero(m: [*c]gsl_matrix_uint) void;
pub extern fn gsl_matrix_uint_set_identity(m: [*c]gsl_matrix_uint) void;
pub extern fn gsl_matrix_uint_set_all(m: [*c]gsl_matrix_uint, x: c_uint) void;
pub extern fn gsl_matrix_uint_fread(stream: ?*FILE, m: [*c]gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_uint, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_uint_memcpy(dest: [*c]gsl_matrix_uint, src: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_swap(m1: [*c]gsl_matrix_uint, m2: [*c]gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_uint, src: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_swap_rows(m: [*c]gsl_matrix_uint, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uint_swap_columns(m: [*c]gsl_matrix_uint, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uint_swap_rowcol(m: [*c]gsl_matrix_uint, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uint_transpose(m: [*c]gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_transpose_memcpy(dest: [*c]gsl_matrix_uint, src: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_uint, src: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_max(m: [*c]const gsl_matrix_uint) c_uint;
pub extern fn gsl_matrix_uint_min(m: [*c]const gsl_matrix_uint) c_uint;
pub extern fn gsl_matrix_uint_minmax(m: [*c]const gsl_matrix_uint, min_out: [*c]c_uint, max_out: [*c]c_uint) void;
pub extern fn gsl_matrix_uint_max_index(m: [*c]const gsl_matrix_uint, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_uint_min_index(m: [*c]const gsl_matrix_uint, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_uint_minmax_index(m: [*c]const gsl_matrix_uint, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_uint_equal(a: [*c]const gsl_matrix_uint, b: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_isnull(m: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_ispos(m: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_isneg(m: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_isnonneg(m: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_norm1(m: [*c]const gsl_matrix_uint) c_uint;
pub extern fn gsl_matrix_uint_add(a: [*c]gsl_matrix_uint, b: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_sub(a: [*c]gsl_matrix_uint, b: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_mul_elements(a: [*c]gsl_matrix_uint, b: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_div_elements(a: [*c]gsl_matrix_uint, b: [*c]const gsl_matrix_uint) c_int;
pub extern fn gsl_matrix_uint_scale(a: [*c]gsl_matrix_uint, x: c_uint) c_int;
pub extern fn gsl_matrix_uint_scale_rows(a: [*c]gsl_matrix_uint, x: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_matrix_uint_scale_columns(a: [*c]gsl_matrix_uint, x: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_matrix_uint_add_constant(a: [*c]gsl_matrix_uint, x: c_uint) c_int;
pub extern fn gsl_matrix_uint_add_diagonal(a: [*c]gsl_matrix_uint, x: c_uint) c_int;
pub extern fn gsl_matrix_uint_get_row(v: [*c]gsl_vector_uint, m: [*c]const gsl_matrix_uint, i: usize) c_int;
pub extern fn gsl_matrix_uint_get_col(v: [*c]gsl_vector_uint, m: [*c]const gsl_matrix_uint, j: usize) c_int;
pub extern fn gsl_matrix_uint_set_row(m: [*c]gsl_matrix_uint, i: usize, v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_matrix_uint_set_col(m: [*c]gsl_matrix_uint, j: usize, v: [*c]const gsl_vector_uint) c_int;
pub extern fn gsl_matrix_uint_get(m: [*c]const gsl_matrix_uint, i: usize, j: usize) c_uint;
pub extern fn gsl_matrix_uint_set(m: [*c]gsl_matrix_uint, i: usize, j: usize, x: c_uint) void;
pub extern fn gsl_matrix_uint_ptr(m: [*c]gsl_matrix_uint, i: usize, j: usize) [*c]c_uint;
pub extern fn gsl_matrix_uint_const_ptr(m: [*c]const gsl_matrix_uint, i: usize, j: usize) [*c]const c_uint;
pub const gsl_matrix_int = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_int = @import("std").mem.zeroes([*c]c_int),
    block: [*c]gsl_block_int = @import("std").mem.zeroes([*c]gsl_block_int),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_int_view = extern struct {
    matrix: gsl_matrix_int = @import("std").mem.zeroes(gsl_matrix_int),
};
pub const gsl_matrix_int_view = _gsl_matrix_int_view;
pub const _gsl_matrix_int_const_view = extern struct {
    matrix: gsl_matrix_int = @import("std").mem.zeroes(gsl_matrix_int),
};
pub const gsl_matrix_int_const_view = _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_alloc(n1: usize, n2: usize) [*c]gsl_matrix_int;
pub extern fn gsl_matrix_int_calloc(n1: usize, n2: usize) [*c]gsl_matrix_int;
pub extern fn gsl_matrix_int_alloc_from_block(b: [*c]gsl_block_int, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_int;
pub extern fn gsl_matrix_int_alloc_from_matrix(m: [*c]gsl_matrix_int, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_int;
pub extern fn gsl_vector_int_alloc_row_from_matrix(m: [*c]gsl_matrix_int, i: usize) [*c]gsl_vector_int;
pub extern fn gsl_vector_int_alloc_col_from_matrix(m: [*c]gsl_matrix_int, j: usize) [*c]gsl_vector_int;
pub extern fn gsl_matrix_int_free(m: [*c]gsl_matrix_int) void;
pub extern fn gsl_matrix_int_submatrix(m: [*c]gsl_matrix_int, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_int_view;
pub extern fn gsl_matrix_int_row(m: [*c]gsl_matrix_int, i: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_column(m: [*c]gsl_matrix_int, j: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_diagonal(m: [*c]gsl_matrix_int) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_subdiagonal(m: [*c]gsl_matrix_int, k: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_superdiagonal(m: [*c]gsl_matrix_int, k: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_subrow(m: [*c]gsl_matrix_int, i: usize, offset: usize, n: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_subcolumn(m: [*c]gsl_matrix_int, j: usize, offset: usize, n: usize) _gsl_vector_int_view;
pub extern fn gsl_matrix_int_view_array(base: [*c]c_int, n1: usize, n2: usize) _gsl_matrix_int_view;
pub extern fn gsl_matrix_int_view_array_with_tda(base: [*c]c_int, n1: usize, n2: usize, tda: usize) _gsl_matrix_int_view;
pub extern fn gsl_matrix_int_view_vector(v: [*c]gsl_vector_int, n1: usize, n2: usize) _gsl_matrix_int_view;
pub extern fn gsl_matrix_int_view_vector_with_tda(v: [*c]gsl_vector_int, n1: usize, n2: usize, tda: usize) _gsl_matrix_int_view;
pub extern fn gsl_matrix_int_const_submatrix(m: [*c]const gsl_matrix_int, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_const_row(m: [*c]const gsl_matrix_int, i: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_column(m: [*c]const gsl_matrix_int, j: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_diagonal(m: [*c]const gsl_matrix_int) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_subdiagonal(m: [*c]const gsl_matrix_int, k: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_superdiagonal(m: [*c]const gsl_matrix_int, k: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_subrow(m: [*c]const gsl_matrix_int, i: usize, offset: usize, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_subcolumn(m: [*c]const gsl_matrix_int, j: usize, offset: usize, n: usize) _gsl_vector_int_const_view;
pub extern fn gsl_matrix_int_const_view_array(base: [*c]const c_int, n1: usize, n2: usize) _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_const_view_array_with_tda(base: [*c]const c_int, n1: usize, n2: usize, tda: usize) _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_const_view_vector(v: [*c]const gsl_vector_int, n1: usize, n2: usize) _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_const_view_vector_with_tda(v: [*c]const gsl_vector_int, n1: usize, n2: usize, tda: usize) _gsl_matrix_int_const_view;
pub extern fn gsl_matrix_int_set_zero(m: [*c]gsl_matrix_int) void;
pub extern fn gsl_matrix_int_set_identity(m: [*c]gsl_matrix_int) void;
pub extern fn gsl_matrix_int_set_all(m: [*c]gsl_matrix_int, x: c_int) void;
pub extern fn gsl_matrix_int_fread(stream: ?*FILE, m: [*c]gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_int, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_int_memcpy(dest: [*c]gsl_matrix_int, src: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_swap(m1: [*c]gsl_matrix_int, m2: [*c]gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_int, src: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_swap_rows(m: [*c]gsl_matrix_int, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_int_swap_columns(m: [*c]gsl_matrix_int, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_int_swap_rowcol(m: [*c]gsl_matrix_int, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_int_transpose(m: [*c]gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_transpose_memcpy(dest: [*c]gsl_matrix_int, src: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_int, src: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_max(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_min(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_minmax(m: [*c]const gsl_matrix_int, min_out: [*c]c_int, max_out: [*c]c_int) void;
pub extern fn gsl_matrix_int_max_index(m: [*c]const gsl_matrix_int, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_int_min_index(m: [*c]const gsl_matrix_int, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_int_minmax_index(m: [*c]const gsl_matrix_int, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_int_equal(a: [*c]const gsl_matrix_int, b: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_isnull(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_ispos(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_isneg(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_isnonneg(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_norm1(m: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_add(a: [*c]gsl_matrix_int, b: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_sub(a: [*c]gsl_matrix_int, b: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_mul_elements(a: [*c]gsl_matrix_int, b: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_div_elements(a: [*c]gsl_matrix_int, b: [*c]const gsl_matrix_int) c_int;
pub extern fn gsl_matrix_int_scale(a: [*c]gsl_matrix_int, x: c_int) c_int;
pub extern fn gsl_matrix_int_scale_rows(a: [*c]gsl_matrix_int, x: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_matrix_int_scale_columns(a: [*c]gsl_matrix_int, x: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_matrix_int_add_constant(a: [*c]gsl_matrix_int, x: c_int) c_int;
pub extern fn gsl_matrix_int_add_diagonal(a: [*c]gsl_matrix_int, x: c_int) c_int;
pub extern fn gsl_matrix_int_get_row(v: [*c]gsl_vector_int, m: [*c]const gsl_matrix_int, i: usize) c_int;
pub extern fn gsl_matrix_int_get_col(v: [*c]gsl_vector_int, m: [*c]const gsl_matrix_int, j: usize) c_int;
pub extern fn gsl_matrix_int_set_row(m: [*c]gsl_matrix_int, i: usize, v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_matrix_int_set_col(m: [*c]gsl_matrix_int, j: usize, v: [*c]const gsl_vector_int) c_int;
pub extern fn gsl_matrix_int_get(m: [*c]const gsl_matrix_int, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_int_set(m: [*c]gsl_matrix_int, i: usize, j: usize, x: c_int) void;
pub extern fn gsl_matrix_int_ptr(m: [*c]gsl_matrix_int, i: usize, j: usize) [*c]c_int;
pub extern fn gsl_matrix_int_const_ptr(m: [*c]const gsl_matrix_int, i: usize, j: usize) [*c]const c_int;
pub const gsl_matrix_ushort = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_ushort = @import("std").mem.zeroes([*c]c_ushort),
    block: [*c]gsl_block_ushort = @import("std").mem.zeroes([*c]gsl_block_ushort),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_ushort_view = extern struct {
    matrix: gsl_matrix_ushort = @import("std").mem.zeroes(gsl_matrix_ushort),
};
pub const gsl_matrix_ushort_view = _gsl_matrix_ushort_view;
pub const _gsl_matrix_ushort_const_view = extern struct {
    matrix: gsl_matrix_ushort = @import("std").mem.zeroes(gsl_matrix_ushort),
};
pub const gsl_matrix_ushort_const_view = _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_alloc(n1: usize, n2: usize) [*c]gsl_matrix_ushort;
pub extern fn gsl_matrix_ushort_calloc(n1: usize, n2: usize) [*c]gsl_matrix_ushort;
pub extern fn gsl_matrix_ushort_alloc_from_block(b: [*c]gsl_block_ushort, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_ushort;
pub extern fn gsl_matrix_ushort_alloc_from_matrix(m: [*c]gsl_matrix_ushort, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_ushort;
pub extern fn gsl_vector_ushort_alloc_row_from_matrix(m: [*c]gsl_matrix_ushort, i: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_vector_ushort_alloc_col_from_matrix(m: [*c]gsl_matrix_ushort, j: usize) [*c]gsl_vector_ushort;
pub extern fn gsl_matrix_ushort_free(m: [*c]gsl_matrix_ushort) void;
pub extern fn gsl_matrix_ushort_submatrix(m: [*c]gsl_matrix_ushort, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_ushort_view;
pub extern fn gsl_matrix_ushort_row(m: [*c]gsl_matrix_ushort, i: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_column(m: [*c]gsl_matrix_ushort, j: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_diagonal(m: [*c]gsl_matrix_ushort) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_subdiagonal(m: [*c]gsl_matrix_ushort, k: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_superdiagonal(m: [*c]gsl_matrix_ushort, k: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_subrow(m: [*c]gsl_matrix_ushort, i: usize, offset: usize, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_subcolumn(m: [*c]gsl_matrix_ushort, j: usize, offset: usize, n: usize) _gsl_vector_ushort_view;
pub extern fn gsl_matrix_ushort_view_array(base: [*c]c_ushort, n1: usize, n2: usize) _gsl_matrix_ushort_view;
pub extern fn gsl_matrix_ushort_view_array_with_tda(base: [*c]c_ushort, n1: usize, n2: usize, tda: usize) _gsl_matrix_ushort_view;
pub extern fn gsl_matrix_ushort_view_vector(v: [*c]gsl_vector_ushort, n1: usize, n2: usize) _gsl_matrix_ushort_view;
pub extern fn gsl_matrix_ushort_view_vector_with_tda(v: [*c]gsl_vector_ushort, n1: usize, n2: usize, tda: usize) _gsl_matrix_ushort_view;
pub extern fn gsl_matrix_ushort_const_submatrix(m: [*c]const gsl_matrix_ushort, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_row(m: [*c]const gsl_matrix_ushort, i: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_column(m: [*c]const gsl_matrix_ushort, j: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_diagonal(m: [*c]const gsl_matrix_ushort) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_subdiagonal(m: [*c]const gsl_matrix_ushort, k: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_superdiagonal(m: [*c]const gsl_matrix_ushort, k: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_subrow(m: [*c]const gsl_matrix_ushort, i: usize, offset: usize, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_subcolumn(m: [*c]const gsl_matrix_ushort, j: usize, offset: usize, n: usize) _gsl_vector_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_view_array(base: [*c]const c_ushort, n1: usize, n2: usize) _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_view_array_with_tda(base: [*c]const c_ushort, n1: usize, n2: usize, tda: usize) _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_view_vector(v: [*c]const gsl_vector_ushort, n1: usize, n2: usize) _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_const_view_vector_with_tda(v: [*c]const gsl_vector_ushort, n1: usize, n2: usize, tda: usize) _gsl_matrix_ushort_const_view;
pub extern fn gsl_matrix_ushort_set_zero(m: [*c]gsl_matrix_ushort) void;
pub extern fn gsl_matrix_ushort_set_identity(m: [*c]gsl_matrix_ushort) void;
pub extern fn gsl_matrix_ushort_set_all(m: [*c]gsl_matrix_ushort, x: c_ushort) void;
pub extern fn gsl_matrix_ushort_fread(stream: ?*FILE, m: [*c]gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_ushort, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_ushort_memcpy(dest: [*c]gsl_matrix_ushort, src: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_swap(m1: [*c]gsl_matrix_ushort, m2: [*c]gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_ushort, src: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_swap_rows(m: [*c]gsl_matrix_ushort, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ushort_swap_columns(m: [*c]gsl_matrix_ushort, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ushort_swap_rowcol(m: [*c]gsl_matrix_ushort, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_ushort_transpose(m: [*c]gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_transpose_memcpy(dest: [*c]gsl_matrix_ushort, src: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_ushort, src: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_max(m: [*c]const gsl_matrix_ushort) c_ushort;
pub extern fn gsl_matrix_ushort_min(m: [*c]const gsl_matrix_ushort) c_ushort;
pub extern fn gsl_matrix_ushort_minmax(m: [*c]const gsl_matrix_ushort, min_out: [*c]c_ushort, max_out: [*c]c_ushort) void;
pub extern fn gsl_matrix_ushort_max_index(m: [*c]const gsl_matrix_ushort, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_ushort_min_index(m: [*c]const gsl_matrix_ushort, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_ushort_minmax_index(m: [*c]const gsl_matrix_ushort, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_ushort_equal(a: [*c]const gsl_matrix_ushort, b: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_isnull(m: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_ispos(m: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_isneg(m: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_isnonneg(m: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_norm1(m: [*c]const gsl_matrix_ushort) c_ushort;
pub extern fn gsl_matrix_ushort_add(a: [*c]gsl_matrix_ushort, b: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_sub(a: [*c]gsl_matrix_ushort, b: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_mul_elements(a: [*c]gsl_matrix_ushort, b: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_div_elements(a: [*c]gsl_matrix_ushort, b: [*c]const gsl_matrix_ushort) c_int;
pub extern fn gsl_matrix_ushort_scale(a: [*c]gsl_matrix_ushort, x: c_ushort) c_int;
pub extern fn gsl_matrix_ushort_scale_rows(a: [*c]gsl_matrix_ushort, x: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_matrix_ushort_scale_columns(a: [*c]gsl_matrix_ushort, x: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_matrix_ushort_add_constant(a: [*c]gsl_matrix_ushort, x: c_ushort) c_int;
pub extern fn gsl_matrix_ushort_add_diagonal(a: [*c]gsl_matrix_ushort, x: c_ushort) c_int;
pub extern fn gsl_matrix_ushort_get_row(v: [*c]gsl_vector_ushort, m: [*c]const gsl_matrix_ushort, i: usize) c_int;
pub extern fn gsl_matrix_ushort_get_col(v: [*c]gsl_vector_ushort, m: [*c]const gsl_matrix_ushort, j: usize) c_int;
pub extern fn gsl_matrix_ushort_set_row(m: [*c]gsl_matrix_ushort, i: usize, v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_matrix_ushort_set_col(m: [*c]gsl_matrix_ushort, j: usize, v: [*c]const gsl_vector_ushort) c_int;
pub extern fn gsl_matrix_ushort_get(m: [*c]const gsl_matrix_ushort, i: usize, j: usize) c_ushort;
pub extern fn gsl_matrix_ushort_set(m: [*c]gsl_matrix_ushort, i: usize, j: usize, x: c_ushort) void;
pub extern fn gsl_matrix_ushort_ptr(m: [*c]gsl_matrix_ushort, i: usize, j: usize) [*c]c_ushort;
pub extern fn gsl_matrix_ushort_const_ptr(m: [*c]const gsl_matrix_ushort, i: usize, j: usize) [*c]const c_ushort;
pub const gsl_matrix_short = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]c_short = @import("std").mem.zeroes([*c]c_short),
    block: [*c]gsl_block_short = @import("std").mem.zeroes([*c]gsl_block_short),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_short_view = extern struct {
    matrix: gsl_matrix_short = @import("std").mem.zeroes(gsl_matrix_short),
};
pub const gsl_matrix_short_view = _gsl_matrix_short_view;
pub const _gsl_matrix_short_const_view = extern struct {
    matrix: gsl_matrix_short = @import("std").mem.zeroes(gsl_matrix_short),
};
pub const gsl_matrix_short_const_view = _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_alloc(n1: usize, n2: usize) [*c]gsl_matrix_short;
pub extern fn gsl_matrix_short_calloc(n1: usize, n2: usize) [*c]gsl_matrix_short;
pub extern fn gsl_matrix_short_alloc_from_block(b: [*c]gsl_block_short, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_short;
pub extern fn gsl_matrix_short_alloc_from_matrix(m: [*c]gsl_matrix_short, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_short;
pub extern fn gsl_vector_short_alloc_row_from_matrix(m: [*c]gsl_matrix_short, i: usize) [*c]gsl_vector_short;
pub extern fn gsl_vector_short_alloc_col_from_matrix(m: [*c]gsl_matrix_short, j: usize) [*c]gsl_vector_short;
pub extern fn gsl_matrix_short_free(m: [*c]gsl_matrix_short) void;
pub extern fn gsl_matrix_short_submatrix(m: [*c]gsl_matrix_short, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_short_view;
pub extern fn gsl_matrix_short_row(m: [*c]gsl_matrix_short, i: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_column(m: [*c]gsl_matrix_short, j: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_diagonal(m: [*c]gsl_matrix_short) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_subdiagonal(m: [*c]gsl_matrix_short, k: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_superdiagonal(m: [*c]gsl_matrix_short, k: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_subrow(m: [*c]gsl_matrix_short, i: usize, offset: usize, n: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_subcolumn(m: [*c]gsl_matrix_short, j: usize, offset: usize, n: usize) _gsl_vector_short_view;
pub extern fn gsl_matrix_short_view_array(base: [*c]c_short, n1: usize, n2: usize) _gsl_matrix_short_view;
pub extern fn gsl_matrix_short_view_array_with_tda(base: [*c]c_short, n1: usize, n2: usize, tda: usize) _gsl_matrix_short_view;
pub extern fn gsl_matrix_short_view_vector(v: [*c]gsl_vector_short, n1: usize, n2: usize) _gsl_matrix_short_view;
pub extern fn gsl_matrix_short_view_vector_with_tda(v: [*c]gsl_vector_short, n1: usize, n2: usize, tda: usize) _gsl_matrix_short_view;
pub extern fn gsl_matrix_short_const_submatrix(m: [*c]const gsl_matrix_short, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_const_row(m: [*c]const gsl_matrix_short, i: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_column(m: [*c]const gsl_matrix_short, j: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_diagonal(m: [*c]const gsl_matrix_short) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_subdiagonal(m: [*c]const gsl_matrix_short, k: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_superdiagonal(m: [*c]const gsl_matrix_short, k: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_subrow(m: [*c]const gsl_matrix_short, i: usize, offset: usize, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_subcolumn(m: [*c]const gsl_matrix_short, j: usize, offset: usize, n: usize) _gsl_vector_short_const_view;
pub extern fn gsl_matrix_short_const_view_array(base: [*c]const c_short, n1: usize, n2: usize) _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_const_view_array_with_tda(base: [*c]const c_short, n1: usize, n2: usize, tda: usize) _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_const_view_vector(v: [*c]const gsl_vector_short, n1: usize, n2: usize) _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_const_view_vector_with_tda(v: [*c]const gsl_vector_short, n1: usize, n2: usize, tda: usize) _gsl_matrix_short_const_view;
pub extern fn gsl_matrix_short_set_zero(m: [*c]gsl_matrix_short) void;
pub extern fn gsl_matrix_short_set_identity(m: [*c]gsl_matrix_short) void;
pub extern fn gsl_matrix_short_set_all(m: [*c]gsl_matrix_short, x: c_short) void;
pub extern fn gsl_matrix_short_fread(stream: ?*FILE, m: [*c]gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_short, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_short_memcpy(dest: [*c]gsl_matrix_short, src: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_swap(m1: [*c]gsl_matrix_short, m2: [*c]gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_short, src: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_swap_rows(m: [*c]gsl_matrix_short, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_short_swap_columns(m: [*c]gsl_matrix_short, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_short_swap_rowcol(m: [*c]gsl_matrix_short, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_short_transpose(m: [*c]gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_transpose_memcpy(dest: [*c]gsl_matrix_short, src: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_short, src: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_max(m: [*c]const gsl_matrix_short) c_short;
pub extern fn gsl_matrix_short_min(m: [*c]const gsl_matrix_short) c_short;
pub extern fn gsl_matrix_short_minmax(m: [*c]const gsl_matrix_short, min_out: [*c]c_short, max_out: [*c]c_short) void;
pub extern fn gsl_matrix_short_max_index(m: [*c]const gsl_matrix_short, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_short_min_index(m: [*c]const gsl_matrix_short, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_short_minmax_index(m: [*c]const gsl_matrix_short, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_short_equal(a: [*c]const gsl_matrix_short, b: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_isnull(m: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_ispos(m: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_isneg(m: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_isnonneg(m: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_norm1(m: [*c]const gsl_matrix_short) c_short;
pub extern fn gsl_matrix_short_add(a: [*c]gsl_matrix_short, b: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_sub(a: [*c]gsl_matrix_short, b: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_mul_elements(a: [*c]gsl_matrix_short, b: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_div_elements(a: [*c]gsl_matrix_short, b: [*c]const gsl_matrix_short) c_int;
pub extern fn gsl_matrix_short_scale(a: [*c]gsl_matrix_short, x: c_short) c_int;
pub extern fn gsl_matrix_short_scale_rows(a: [*c]gsl_matrix_short, x: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_matrix_short_scale_columns(a: [*c]gsl_matrix_short, x: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_matrix_short_add_constant(a: [*c]gsl_matrix_short, x: c_short) c_int;
pub extern fn gsl_matrix_short_add_diagonal(a: [*c]gsl_matrix_short, x: c_short) c_int;
pub extern fn gsl_matrix_short_get_row(v: [*c]gsl_vector_short, m: [*c]const gsl_matrix_short, i: usize) c_int;
pub extern fn gsl_matrix_short_get_col(v: [*c]gsl_vector_short, m: [*c]const gsl_matrix_short, j: usize) c_int;
pub extern fn gsl_matrix_short_set_row(m: [*c]gsl_matrix_short, i: usize, v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_matrix_short_set_col(m: [*c]gsl_matrix_short, j: usize, v: [*c]const gsl_vector_short) c_int;
pub extern fn gsl_matrix_short_get(m: [*c]const gsl_matrix_short, i: usize, j: usize) c_short;
pub extern fn gsl_matrix_short_set(m: [*c]gsl_matrix_short, i: usize, j: usize, x: c_short) void;
pub extern fn gsl_matrix_short_ptr(m: [*c]gsl_matrix_short, i: usize, j: usize) [*c]c_short;
pub extern fn gsl_matrix_short_const_ptr(m: [*c]const gsl_matrix_short, i: usize, j: usize) [*c]const c_short;
pub const gsl_matrix_uchar = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
    block: [*c]gsl_block_uchar = @import("std").mem.zeroes([*c]gsl_block_uchar),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_uchar_view = extern struct {
    matrix: gsl_matrix_uchar = @import("std").mem.zeroes(gsl_matrix_uchar),
};
pub const gsl_matrix_uchar_view = _gsl_matrix_uchar_view;
pub const _gsl_matrix_uchar_const_view = extern struct {
    matrix: gsl_matrix_uchar = @import("std").mem.zeroes(gsl_matrix_uchar),
};
pub const gsl_matrix_uchar_const_view = _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_alloc(n1: usize, n2: usize) [*c]gsl_matrix_uchar;
pub extern fn gsl_matrix_uchar_calloc(n1: usize, n2: usize) [*c]gsl_matrix_uchar;
pub extern fn gsl_matrix_uchar_alloc_from_block(b: [*c]gsl_block_uchar, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_uchar;
pub extern fn gsl_matrix_uchar_alloc_from_matrix(m: [*c]gsl_matrix_uchar, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_uchar;
pub extern fn gsl_vector_uchar_alloc_row_from_matrix(m: [*c]gsl_matrix_uchar, i: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_vector_uchar_alloc_col_from_matrix(m: [*c]gsl_matrix_uchar, j: usize) [*c]gsl_vector_uchar;
pub extern fn gsl_matrix_uchar_free(m: [*c]gsl_matrix_uchar) void;
pub extern fn gsl_matrix_uchar_submatrix(m: [*c]gsl_matrix_uchar, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_uchar_view;
pub extern fn gsl_matrix_uchar_row(m: [*c]gsl_matrix_uchar, i: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_column(m: [*c]gsl_matrix_uchar, j: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_diagonal(m: [*c]gsl_matrix_uchar) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_subdiagonal(m: [*c]gsl_matrix_uchar, k: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_superdiagonal(m: [*c]gsl_matrix_uchar, k: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_subrow(m: [*c]gsl_matrix_uchar, i: usize, offset: usize, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_subcolumn(m: [*c]gsl_matrix_uchar, j: usize, offset: usize, n: usize) _gsl_vector_uchar_view;
pub extern fn gsl_matrix_uchar_view_array(base: [*c]u8, n1: usize, n2: usize) _gsl_matrix_uchar_view;
pub extern fn gsl_matrix_uchar_view_array_with_tda(base: [*c]u8, n1: usize, n2: usize, tda: usize) _gsl_matrix_uchar_view;
pub extern fn gsl_matrix_uchar_view_vector(v: [*c]gsl_vector_uchar, n1: usize, n2: usize) _gsl_matrix_uchar_view;
pub extern fn gsl_matrix_uchar_view_vector_with_tda(v: [*c]gsl_vector_uchar, n1: usize, n2: usize, tda: usize) _gsl_matrix_uchar_view;
pub extern fn gsl_matrix_uchar_const_submatrix(m: [*c]const gsl_matrix_uchar, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_row(m: [*c]const gsl_matrix_uchar, i: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_column(m: [*c]const gsl_matrix_uchar, j: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_diagonal(m: [*c]const gsl_matrix_uchar) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_subdiagonal(m: [*c]const gsl_matrix_uchar, k: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_superdiagonal(m: [*c]const gsl_matrix_uchar, k: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_subrow(m: [*c]const gsl_matrix_uchar, i: usize, offset: usize, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_subcolumn(m: [*c]const gsl_matrix_uchar, j: usize, offset: usize, n: usize) _gsl_vector_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_view_array(base: [*c]const u8, n1: usize, n2: usize) _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_view_array_with_tda(base: [*c]const u8, n1: usize, n2: usize, tda: usize) _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_view_vector(v: [*c]const gsl_vector_uchar, n1: usize, n2: usize) _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_const_view_vector_with_tda(v: [*c]const gsl_vector_uchar, n1: usize, n2: usize, tda: usize) _gsl_matrix_uchar_const_view;
pub extern fn gsl_matrix_uchar_set_zero(m: [*c]gsl_matrix_uchar) void;
pub extern fn gsl_matrix_uchar_set_identity(m: [*c]gsl_matrix_uchar) void;
pub extern fn gsl_matrix_uchar_set_all(m: [*c]gsl_matrix_uchar, x: u8) void;
pub extern fn gsl_matrix_uchar_fread(stream: ?*FILE, m: [*c]gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_uchar, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_uchar_memcpy(dest: [*c]gsl_matrix_uchar, src: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_swap(m1: [*c]gsl_matrix_uchar, m2: [*c]gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_uchar, src: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_swap_rows(m: [*c]gsl_matrix_uchar, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uchar_swap_columns(m: [*c]gsl_matrix_uchar, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uchar_swap_rowcol(m: [*c]gsl_matrix_uchar, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_uchar_transpose(m: [*c]gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_transpose_memcpy(dest: [*c]gsl_matrix_uchar, src: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_uchar, src: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_max(m: [*c]const gsl_matrix_uchar) u8;
pub extern fn gsl_matrix_uchar_min(m: [*c]const gsl_matrix_uchar) u8;
pub extern fn gsl_matrix_uchar_minmax(m: [*c]const gsl_matrix_uchar, min_out: [*c]u8, max_out: [*c]u8) void;
pub extern fn gsl_matrix_uchar_max_index(m: [*c]const gsl_matrix_uchar, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_uchar_min_index(m: [*c]const gsl_matrix_uchar, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_uchar_minmax_index(m: [*c]const gsl_matrix_uchar, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_uchar_equal(a: [*c]const gsl_matrix_uchar, b: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_isnull(m: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_ispos(m: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_isneg(m: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_isnonneg(m: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_norm1(m: [*c]const gsl_matrix_uchar) u8;
pub extern fn gsl_matrix_uchar_add(a: [*c]gsl_matrix_uchar, b: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_sub(a: [*c]gsl_matrix_uchar, b: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_mul_elements(a: [*c]gsl_matrix_uchar, b: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_div_elements(a: [*c]gsl_matrix_uchar, b: [*c]const gsl_matrix_uchar) c_int;
pub extern fn gsl_matrix_uchar_scale(a: [*c]gsl_matrix_uchar, x: u8) c_int;
pub extern fn gsl_matrix_uchar_scale_rows(a: [*c]gsl_matrix_uchar, x: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_matrix_uchar_scale_columns(a: [*c]gsl_matrix_uchar, x: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_matrix_uchar_add_constant(a: [*c]gsl_matrix_uchar, x: u8) c_int;
pub extern fn gsl_matrix_uchar_add_diagonal(a: [*c]gsl_matrix_uchar, x: u8) c_int;
pub extern fn gsl_matrix_uchar_get_row(v: [*c]gsl_vector_uchar, m: [*c]const gsl_matrix_uchar, i: usize) c_int;
pub extern fn gsl_matrix_uchar_get_col(v: [*c]gsl_vector_uchar, m: [*c]const gsl_matrix_uchar, j: usize) c_int;
pub extern fn gsl_matrix_uchar_set_row(m: [*c]gsl_matrix_uchar, i: usize, v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_matrix_uchar_set_col(m: [*c]gsl_matrix_uchar, j: usize, v: [*c]const gsl_vector_uchar) c_int;
pub extern fn gsl_matrix_uchar_get(m: [*c]const gsl_matrix_uchar, i: usize, j: usize) u8;
pub extern fn gsl_matrix_uchar_set(m: [*c]gsl_matrix_uchar, i: usize, j: usize, x: u8) void;
pub extern fn gsl_matrix_uchar_ptr(m: [*c]gsl_matrix_uchar, i: usize, j: usize) [*c]u8;
pub extern fn gsl_matrix_uchar_const_ptr(m: [*c]const gsl_matrix_uchar, i: usize, j: usize) [*c]const u8;
pub const gsl_matrix_char = extern struct {
    size1: usize = @import("std").mem.zeroes(usize),
    size2: usize = @import("std").mem.zeroes(usize),
    tda: usize = @import("std").mem.zeroes(usize),
    data: [*c]u8 = @import("std").mem.zeroes([*c]u8),
    block: [*c]gsl_block_char = @import("std").mem.zeroes([*c]gsl_block_char),
    owner: c_int = @import("std").mem.zeroes(c_int),
};
pub const _gsl_matrix_char_view = extern struct {
    matrix: gsl_matrix_char = @import("std").mem.zeroes(gsl_matrix_char),
};
pub const gsl_matrix_char_view = _gsl_matrix_char_view;
pub const _gsl_matrix_char_const_view = extern struct {
    matrix: gsl_matrix_char = @import("std").mem.zeroes(gsl_matrix_char),
};
pub const gsl_matrix_char_const_view = _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_alloc(n1: usize, n2: usize) [*c]gsl_matrix_char;
pub extern fn gsl_matrix_char_calloc(n1: usize, n2: usize) [*c]gsl_matrix_char;
pub extern fn gsl_matrix_char_alloc_from_block(b: [*c]gsl_block_char, offset: usize, n1: usize, n2: usize, d2: usize) [*c]gsl_matrix_char;
pub extern fn gsl_matrix_char_alloc_from_matrix(m: [*c]gsl_matrix_char, k1: usize, k2: usize, n1: usize, n2: usize) [*c]gsl_matrix_char;
pub extern fn gsl_vector_char_alloc_row_from_matrix(m: [*c]gsl_matrix_char, i: usize) [*c]gsl_vector_char;
pub extern fn gsl_vector_char_alloc_col_from_matrix(m: [*c]gsl_matrix_char, j: usize) [*c]gsl_vector_char;
pub extern fn gsl_matrix_char_free(m: [*c]gsl_matrix_char) void;
pub extern fn gsl_matrix_char_submatrix(m: [*c]gsl_matrix_char, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_char_view;
pub extern fn gsl_matrix_char_row(m: [*c]gsl_matrix_char, i: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_column(m: [*c]gsl_matrix_char, j: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_diagonal(m: [*c]gsl_matrix_char) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_subdiagonal(m: [*c]gsl_matrix_char, k: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_superdiagonal(m: [*c]gsl_matrix_char, k: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_subrow(m: [*c]gsl_matrix_char, i: usize, offset: usize, n: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_subcolumn(m: [*c]gsl_matrix_char, j: usize, offset: usize, n: usize) _gsl_vector_char_view;
pub extern fn gsl_matrix_char_view_array(base: [*c]u8, n1: usize, n2: usize) _gsl_matrix_char_view;
pub extern fn gsl_matrix_char_view_array_with_tda(base: [*c]u8, n1: usize, n2: usize, tda: usize) _gsl_matrix_char_view;
pub extern fn gsl_matrix_char_view_vector(v: [*c]gsl_vector_char, n1: usize, n2: usize) _gsl_matrix_char_view;
pub extern fn gsl_matrix_char_view_vector_with_tda(v: [*c]gsl_vector_char, n1: usize, n2: usize, tda: usize) _gsl_matrix_char_view;
pub extern fn gsl_matrix_char_const_submatrix(m: [*c]const gsl_matrix_char, i: usize, j: usize, n1: usize, n2: usize) _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_const_row(m: [*c]const gsl_matrix_char, i: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_column(m: [*c]const gsl_matrix_char, j: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_diagonal(m: [*c]const gsl_matrix_char) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_subdiagonal(m: [*c]const gsl_matrix_char, k: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_superdiagonal(m: [*c]const gsl_matrix_char, k: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_subrow(m: [*c]const gsl_matrix_char, i: usize, offset: usize, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_subcolumn(m: [*c]const gsl_matrix_char, j: usize, offset: usize, n: usize) _gsl_vector_char_const_view;
pub extern fn gsl_matrix_char_const_view_array(base: [*c]const u8, n1: usize, n2: usize) _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_const_view_array_with_tda(base: [*c]const u8, n1: usize, n2: usize, tda: usize) _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_const_view_vector(v: [*c]const gsl_vector_char, n1: usize, n2: usize) _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_const_view_vector_with_tda(v: [*c]const gsl_vector_char, n1: usize, n2: usize, tda: usize) _gsl_matrix_char_const_view;
pub extern fn gsl_matrix_char_set_zero(m: [*c]gsl_matrix_char) void;
pub extern fn gsl_matrix_char_set_identity(m: [*c]gsl_matrix_char) void;
pub extern fn gsl_matrix_char_set_all(m: [*c]gsl_matrix_char, x: u8) void;
pub extern fn gsl_matrix_char_fread(stream: ?*FILE, m: [*c]gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_fwrite(stream: ?*FILE, m: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_fscanf(stream: ?*FILE, m: [*c]gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_fprintf(stream: ?*FILE, m: [*c]const gsl_matrix_char, format: [*c]const u8) c_int;
pub extern fn gsl_matrix_char_memcpy(dest: [*c]gsl_matrix_char, src: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_swap(m1: [*c]gsl_matrix_char, m2: [*c]gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_tricpy(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_char, src: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_swap_rows(m: [*c]gsl_matrix_char, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_char_swap_columns(m: [*c]gsl_matrix_char, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_char_swap_rowcol(m: [*c]gsl_matrix_char, i: usize, j: usize) c_int;
pub extern fn gsl_matrix_char_transpose(m: [*c]gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_transpose_memcpy(dest: [*c]gsl_matrix_char, src: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_transpose_tricpy(Uplo_src: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, dest: [*c]gsl_matrix_char, src: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_max(m: [*c]const gsl_matrix_char) u8;
pub extern fn gsl_matrix_char_min(m: [*c]const gsl_matrix_char) u8;
pub extern fn gsl_matrix_char_minmax(m: [*c]const gsl_matrix_char, min_out: [*c]u8, max_out: [*c]u8) void;
pub extern fn gsl_matrix_char_max_index(m: [*c]const gsl_matrix_char, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_char_min_index(m: [*c]const gsl_matrix_char, imin: [*c]usize, jmin: [*c]usize) void;
pub extern fn gsl_matrix_char_minmax_index(m: [*c]const gsl_matrix_char, imin: [*c]usize, jmin: [*c]usize, imax: [*c]usize, jmax: [*c]usize) void;
pub extern fn gsl_matrix_char_equal(a: [*c]const gsl_matrix_char, b: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_isnull(m: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_ispos(m: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_isneg(m: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_isnonneg(m: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_norm1(m: [*c]const gsl_matrix_char) u8;
pub extern fn gsl_matrix_char_add(a: [*c]gsl_matrix_char, b: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_sub(a: [*c]gsl_matrix_char, b: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_mul_elements(a: [*c]gsl_matrix_char, b: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_div_elements(a: [*c]gsl_matrix_char, b: [*c]const gsl_matrix_char) c_int;
pub extern fn gsl_matrix_char_scale(a: [*c]gsl_matrix_char, x: u8) c_int;
pub extern fn gsl_matrix_char_scale_rows(a: [*c]gsl_matrix_char, x: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_matrix_char_scale_columns(a: [*c]gsl_matrix_char, x: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_matrix_char_add_constant(a: [*c]gsl_matrix_char, x: u8) c_int;
pub extern fn gsl_matrix_char_add_diagonal(a: [*c]gsl_matrix_char, x: u8) c_int;
pub extern fn gsl_matrix_char_get_row(v: [*c]gsl_vector_char, m: [*c]const gsl_matrix_char, i: usize) c_int;
pub extern fn gsl_matrix_char_get_col(v: [*c]gsl_vector_char, m: [*c]const gsl_matrix_char, j: usize) c_int;
pub extern fn gsl_matrix_char_set_row(m: [*c]gsl_matrix_char, i: usize, v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_matrix_char_set_col(m: [*c]gsl_matrix_char, j: usize, v: [*c]const gsl_vector_char) c_int;
pub extern fn gsl_matrix_char_get(m: [*c]const gsl_matrix_char, i: usize, j: usize) u8;
pub extern fn gsl_matrix_char_set(m: [*c]gsl_matrix_char, i: usize, j: usize, x: u8) void;
pub extern fn gsl_matrix_char_ptr(m: [*c]gsl_matrix_char, i: usize, j: usize) [*c]u8;
pub extern fn gsl_matrix_char_const_ptr(m: [*c]const gsl_matrix_char, i: usize, j: usize) [*c]const u8;
pub const float_t = f32;
pub const double_t = f64;
pub extern fn __fpclassify(__value: f64) c_int;
pub extern fn __signbit(__value: f64) c_int;
pub extern fn __isinf(__value: f64) c_int;
pub extern fn __finite(__value: f64) c_int;
pub extern fn __isnan(__value: f64) c_int;
pub extern fn __iseqsig(__x: f64, __y: f64) c_int;
pub extern fn __issignaling(__value: f64) c_int;
pub extern fn acos(__x: f64) f64;
pub extern fn __acos(__x: f64) f64;
pub extern fn asin(__x: f64) f64;
pub extern fn __asin(__x: f64) f64;
pub extern fn atan(__x: f64) f64;
pub extern fn __atan(__x: f64) f64;
pub extern fn atan2(__y: f64, __x: f64) f64;
pub extern fn __atan2(__y: f64, __x: f64) f64;
pub extern fn cos(__x: f64) f64;
pub extern fn __cos(__x: f64) f64;
pub extern fn sin(__x: f64) f64;
pub extern fn __sin(__x: f64) f64;
pub extern fn tan(__x: f64) f64;
pub extern fn __tan(__x: f64) f64;
pub extern fn cosh(__x: f64) f64;
pub extern fn __cosh(__x: f64) f64;
pub extern fn sinh(__x: f64) f64;
pub extern fn __sinh(__x: f64) f64;
pub extern fn tanh(__x: f64) f64;
pub extern fn __tanh(__x: f64) f64;
pub extern fn acosh(__x: f64) f64;
pub extern fn __acosh(__x: f64) f64;
pub extern fn asinh(__x: f64) f64;
pub extern fn __asinh(__x: f64) f64;
pub extern fn atanh(__x: f64) f64;
pub extern fn __atanh(__x: f64) f64;
pub extern fn exp(__x: f64) f64;
pub extern fn __exp(__x: f64) f64;
pub extern fn frexp(__x: f64, __exponent: [*c]c_int) f64;
pub extern fn __frexp(__x: f64, __exponent: [*c]c_int) f64;
pub extern fn ldexp(__x: f64, __exponent: c_int) f64;
pub extern fn __ldexp(__x: f64, __exponent: c_int) f64;
pub extern fn log(__x: f64) f64;
pub extern fn __log(__x: f64) f64;
pub extern fn log10(__x: f64) f64;
pub extern fn __log10(__x: f64) f64;
pub extern fn modf(__x: f64, __iptr: [*c]f64) f64;
pub extern fn __modf(__x: f64, __iptr: [*c]f64) f64;
pub extern fn expm1(__x: f64) f64;
pub extern fn __expm1(__x: f64) f64;
pub extern fn log1p(__x: f64) f64;
pub extern fn __log1p(__x: f64) f64;
pub extern fn logb(__x: f64) f64;
pub extern fn __logb(__x: f64) f64;
pub extern fn exp2(__x: f64) f64;
pub extern fn __exp2(__x: f64) f64;
pub extern fn log2(__x: f64) f64;
pub extern fn __log2(__x: f64) f64;
pub extern fn pow(__x: f64, __y: f64) f64;
pub extern fn __pow(__x: f64, __y: f64) f64;
pub extern fn sqrt(__x: f64) f64;
pub extern fn __sqrt(__x: f64) f64;
pub extern fn hypot(__x: f64, __y: f64) f64;
pub extern fn __hypot(__x: f64, __y: f64) f64;
pub extern fn cbrt(__x: f64) f64;
pub extern fn __cbrt(__x: f64) f64;
pub extern fn ceil(__x: f64) f64;
pub extern fn fabs(__x: f64) f64;
pub extern fn floor(__x: f64) f64;
pub extern fn fmod(__x: f64, __y: f64) f64;
pub extern fn __fmod(__x: f64, __y: f64) f64;
pub extern fn isinf(__value: f64) c_int;
pub extern fn finite(__value: f64) c_int;
pub extern fn drem(__x: f64, __y: f64) f64;
pub extern fn __drem(__x: f64, __y: f64) f64;
pub extern fn significand(__x: f64) f64;
pub extern fn __significand(__x: f64) f64;
pub extern fn copysign(__x: f64, __y: f64) f64;
pub extern fn nan(__tagb: [*c]const u8) f64;
pub extern fn __nan(__tagb: [*c]const u8) f64;
pub extern fn isnan(__value: f64) c_int;
pub extern fn j0(f64) f64;
pub extern fn __j0(f64) f64;
pub extern fn j1(f64) f64;
pub extern fn __j1(f64) f64;
pub extern fn jn(c_int, f64) f64;
pub extern fn __jn(c_int, f64) f64;
pub extern fn y0(f64) f64;
pub extern fn __y0(f64) f64;
pub extern fn y1(f64) f64;
pub extern fn __y1(f64) f64;
pub extern fn yn(c_int, f64) f64;
pub extern fn __yn(c_int, f64) f64;
pub extern fn erf(f64) f64;
pub extern fn __erf(f64) f64;
pub extern fn erfc(f64) f64;
pub extern fn __erfc(f64) f64;
pub extern fn lgamma(f64) f64;
pub extern fn __lgamma(f64) f64;
pub extern fn tgamma(f64) f64;
pub extern fn __tgamma(f64) f64;
pub extern fn gamma(f64) f64;
pub extern fn __gamma(f64) f64;
pub extern fn lgamma_r(f64, __signgamp: [*c]c_int) f64;
pub extern fn __lgamma_r(f64, __signgamp: [*c]c_int) f64;
pub extern fn rint(__x: f64) f64;
pub extern fn __rint(__x: f64) f64;
pub extern fn nextafter(__x: f64, __y: f64) f64;
pub extern fn __nextafter(__x: f64, __y: f64) f64;
pub extern fn nexttoward(__x: f64, __y: c_longdouble) f64;
pub extern fn __nexttoward(__x: f64, __y: c_longdouble) f64;
pub extern fn remainder(__x: f64, __y: f64) f64;
pub extern fn __remainder(__x: f64, __y: f64) f64;
pub extern fn scalbn(__x: f64, __n: c_int) f64;
pub extern fn __scalbn(__x: f64, __n: c_int) f64;
pub extern fn ilogb(__x: f64) c_int;
pub extern fn __ilogb(__x: f64) c_int;
pub extern fn scalbln(__x: f64, __n: c_long) f64;
pub extern fn __scalbln(__x: f64, __n: c_long) f64;
pub extern fn nearbyint(__x: f64) f64;
pub extern fn __nearbyint(__x: f64) f64;
pub extern fn round(__x: f64) f64;
pub extern fn trunc(__x: f64) f64;
pub extern fn remquo(__x: f64, __y: f64, __quo: [*c]c_int) f64;
pub extern fn __remquo(__x: f64, __y: f64, __quo: [*c]c_int) f64;
pub extern fn lrint(__x: f64) c_long;
pub extern fn __lrint(__x: f64) c_long;
pub extern fn llrint(__x: f64) c_longlong;
pub extern fn __llrint(__x: f64) c_longlong;
pub extern fn lround(__x: f64) c_long;
pub extern fn __lround(__x: f64) c_long;
pub extern fn llround(__x: f64) c_longlong;
pub extern fn __llround(__x: f64) c_longlong;
pub extern fn fdim(__x: f64, __y: f64) f64;
pub extern fn __fdim(__x: f64, __y: f64) f64;
pub extern fn fmax(__x: f64, __y: f64) f64;
pub extern fn fmin(__x: f64, __y: f64) f64;
pub extern fn fma(__x: f64, __y: f64, __z: f64) f64;
pub extern fn __fma(__x: f64, __y: f64, __z: f64) f64;
pub extern fn scalb(__x: f64, __n: f64) f64;
pub extern fn __scalb(__x: f64, __n: f64) f64;
pub extern fn __fpclassifyf(__value: f32) c_int;
pub extern fn __signbitf(__value: f32) c_int;
pub extern fn __isinff(__value: f32) c_int;
pub extern fn __finitef(__value: f32) c_int;
pub extern fn __isnanf(__value: f32) c_int;
pub extern fn __iseqsigf(__x: f32, __y: f32) c_int;
pub extern fn __issignalingf(__value: f32) c_int;
pub extern fn acosf(__x: f32) f32;
pub extern fn __acosf(__x: f32) f32;
pub extern fn asinf(__x: f32) f32;
pub extern fn __asinf(__x: f32) f32;
pub extern fn atanf(__x: f32) f32;
pub extern fn __atanf(__x: f32) f32;
pub extern fn atan2f(__y: f32, __x: f32) f32;
pub extern fn __atan2f(__y: f32, __x: f32) f32;
pub extern fn cosf(__x: f32) f32;
pub extern fn __cosf(__x: f32) f32;
pub extern fn sinf(__x: f32) f32;
pub extern fn __sinf(__x: f32) f32;
pub extern fn tanf(__x: f32) f32;
pub extern fn __tanf(__x: f32) f32;
pub extern fn coshf(__x: f32) f32;
pub extern fn __coshf(__x: f32) f32;
pub extern fn sinhf(__x: f32) f32;
pub extern fn __sinhf(__x: f32) f32;
pub extern fn tanhf(__x: f32) f32;
pub extern fn __tanhf(__x: f32) f32;
pub extern fn acoshf(__x: f32) f32;
pub extern fn __acoshf(__x: f32) f32;
pub extern fn asinhf(__x: f32) f32;
pub extern fn __asinhf(__x: f32) f32;
pub extern fn atanhf(__x: f32) f32;
pub extern fn __atanhf(__x: f32) f32;
pub extern fn expf(__x: f32) f32;
pub extern fn __expf(__x: f32) f32;
pub extern fn frexpf(__x: f32, __exponent: [*c]c_int) f32;
pub extern fn __frexpf(__x: f32, __exponent: [*c]c_int) f32;
pub extern fn ldexpf(__x: f32, __exponent: c_int) f32;
pub extern fn __ldexpf(__x: f32, __exponent: c_int) f32;
pub extern fn logf(__x: f32) f32;
pub extern fn __logf(__x: f32) f32;
pub extern fn log10f(__x: f32) f32;
pub extern fn __log10f(__x: f32) f32;
pub extern fn modff(__x: f32, __iptr: [*c]f32) f32;
pub extern fn __modff(__x: f32, __iptr: [*c]f32) f32;
pub extern fn expm1f(__x: f32) f32;
pub extern fn __expm1f(__x: f32) f32;
pub extern fn log1pf(__x: f32) f32;
pub extern fn __log1pf(__x: f32) f32;
pub extern fn logbf(__x: f32) f32;
pub extern fn __logbf(__x: f32) f32;
pub extern fn exp2f(__x: f32) f32;
pub extern fn __exp2f(__x: f32) f32;
pub extern fn log2f(__x: f32) f32;
pub extern fn __log2f(__x: f32) f32;
pub extern fn powf(__x: f32, __y: f32) f32;
pub extern fn __powf(__x: f32, __y: f32) f32;
pub extern fn sqrtf(__x: f32) f32;
pub extern fn __sqrtf(__x: f32) f32;
pub extern fn hypotf(__x: f32, __y: f32) f32;
pub extern fn __hypotf(__x: f32, __y: f32) f32;
pub extern fn cbrtf(__x: f32) f32;
pub extern fn __cbrtf(__x: f32) f32;
pub extern fn ceilf(__x: f32) f32;
pub extern fn fabsf(__x: f32) f32;
pub extern fn floorf(__x: f32) f32;
pub extern fn fmodf(__x: f32, __y: f32) f32;
pub extern fn __fmodf(__x: f32, __y: f32) f32;
pub extern fn isinff(__value: f32) c_int;
pub extern fn finitef(__value: f32) c_int;
pub extern fn dremf(__x: f32, __y: f32) f32;
pub extern fn __dremf(__x: f32, __y: f32) f32;
pub extern fn significandf(__x: f32) f32;
pub extern fn __significandf(__x: f32) f32;
pub extern fn copysignf(__x: f32, __y: f32) f32;
pub extern fn nanf(__tagb: [*c]const u8) f32;
pub extern fn __nanf(__tagb: [*c]const u8) f32;
pub extern fn isnanf(__value: f32) c_int;
pub extern fn j0f(f32) f32;
pub extern fn __j0f(f32) f32;
pub extern fn j1f(f32) f32;
pub extern fn __j1f(f32) f32;
pub extern fn jnf(c_int, f32) f32;
pub extern fn __jnf(c_int, f32) f32;
pub extern fn y0f(f32) f32;
pub extern fn __y0f(f32) f32;
pub extern fn y1f(f32) f32;
pub extern fn __y1f(f32) f32;
pub extern fn ynf(c_int, f32) f32;
pub extern fn __ynf(c_int, f32) f32;
pub extern fn erff(f32) f32;
pub extern fn __erff(f32) f32;
pub extern fn erfcf(f32) f32;
pub extern fn __erfcf(f32) f32;
pub extern fn lgammaf(f32) f32;
pub extern fn __lgammaf(f32) f32;
pub extern fn tgammaf(f32) f32;
pub extern fn __tgammaf(f32) f32;
pub extern fn gammaf(f32) f32;
pub extern fn __gammaf(f32) f32;
pub extern fn lgammaf_r(f32, __signgamp: [*c]c_int) f32;
pub extern fn __lgammaf_r(f32, __signgamp: [*c]c_int) f32;
pub extern fn rintf(__x: f32) f32;
pub extern fn __rintf(__x: f32) f32;
pub extern fn nextafterf(__x: f32, __y: f32) f32;
pub extern fn __nextafterf(__x: f32, __y: f32) f32;
pub extern fn nexttowardf(__x: f32, __y: c_longdouble) f32;
pub extern fn __nexttowardf(__x: f32, __y: c_longdouble) f32;
pub extern fn remainderf(__x: f32, __y: f32) f32;
pub extern fn __remainderf(__x: f32, __y: f32) f32;
pub extern fn scalbnf(__x: f32, __n: c_int) f32;
pub extern fn __scalbnf(__x: f32, __n: c_int) f32;
pub extern fn ilogbf(__x: f32) c_int;
pub extern fn __ilogbf(__x: f32) c_int;
pub extern fn scalblnf(__x: f32, __n: c_long) f32;
pub extern fn __scalblnf(__x: f32, __n: c_long) f32;
pub extern fn nearbyintf(__x: f32) f32;
pub extern fn __nearbyintf(__x: f32) f32;
pub extern fn roundf(__x: f32) f32;
pub extern fn truncf(__x: f32) f32;
pub extern fn remquof(__x: f32, __y: f32, __quo: [*c]c_int) f32;
pub extern fn __remquof(__x: f32, __y: f32, __quo: [*c]c_int) f32;
pub extern fn lrintf(__x: f32) c_long;
pub extern fn __lrintf(__x: f32) c_long;
pub extern fn llrintf(__x: f32) c_longlong;
pub extern fn __llrintf(__x: f32) c_longlong;
pub extern fn lroundf(__x: f32) c_long;
pub extern fn __lroundf(__x: f32) c_long;
pub extern fn llroundf(__x: f32) c_longlong;
pub extern fn __llroundf(__x: f32) c_longlong;
pub extern fn fdimf(__x: f32, __y: f32) f32;
pub extern fn __fdimf(__x: f32, __y: f32) f32;
pub extern fn fmaxf(__x: f32, __y: f32) f32;
pub extern fn fminf(__x: f32, __y: f32) f32;
pub extern fn fmaf(__x: f32, __y: f32, __z: f32) f32;
pub extern fn __fmaf(__x: f32, __y: f32, __z: f32) f32;
pub extern fn scalbf(__x: f32, __n: f32) f32;
pub extern fn __scalbf(__x: f32, __n: f32) f32;
pub extern fn __fpclassifyl(__value: c_longdouble) c_int;
pub extern fn __signbitl(__value: c_longdouble) c_int;
pub extern fn __isinfl(__value: c_longdouble) c_int;
pub extern fn __finitel(__value: c_longdouble) c_int;
pub extern fn __isnanl(__value: c_longdouble) c_int;
pub extern fn __iseqsigl(__x: c_longdouble, __y: c_longdouble) c_int;
pub extern fn __issignalingl(__value: c_longdouble) c_int;
pub extern fn acosl(__x: c_longdouble) c_longdouble;
pub extern fn __acosl(__x: c_longdouble) c_longdouble;
pub extern fn asinl(__x: c_longdouble) c_longdouble;
pub extern fn __asinl(__x: c_longdouble) c_longdouble;
pub extern fn atanl(__x: c_longdouble) c_longdouble;
pub extern fn __atanl(__x: c_longdouble) c_longdouble;
pub extern fn atan2l(__y: c_longdouble, __x: c_longdouble) c_longdouble;
pub extern fn __atan2l(__y: c_longdouble, __x: c_longdouble) c_longdouble;
pub extern fn cosl(__x: c_longdouble) c_longdouble;
pub extern fn __cosl(__x: c_longdouble) c_longdouble;
pub extern fn sinl(__x: c_longdouble) c_longdouble;
pub extern fn __sinl(__x: c_longdouble) c_longdouble;
pub extern fn tanl(__x: c_longdouble) c_longdouble;
pub extern fn __tanl(__x: c_longdouble) c_longdouble;
pub extern fn coshl(__x: c_longdouble) c_longdouble;
pub extern fn __coshl(__x: c_longdouble) c_longdouble;
pub extern fn sinhl(__x: c_longdouble) c_longdouble;
pub extern fn __sinhl(__x: c_longdouble) c_longdouble;
pub extern fn tanhl(__x: c_longdouble) c_longdouble;
pub extern fn __tanhl(__x: c_longdouble) c_longdouble;
pub extern fn acoshl(__x: c_longdouble) c_longdouble;
pub extern fn __acoshl(__x: c_longdouble) c_longdouble;
pub extern fn asinhl(__x: c_longdouble) c_longdouble;
pub extern fn __asinhl(__x: c_longdouble) c_longdouble;
pub extern fn atanhl(__x: c_longdouble) c_longdouble;
pub extern fn __atanhl(__x: c_longdouble) c_longdouble;
pub extern fn expl(__x: c_longdouble) c_longdouble;
pub extern fn __expl(__x: c_longdouble) c_longdouble;
pub extern fn frexpl(__x: c_longdouble, __exponent: [*c]c_int) c_longdouble;
pub extern fn __frexpl(__x: c_longdouble, __exponent: [*c]c_int) c_longdouble;
pub extern fn ldexpl(__x: c_longdouble, __exponent: c_int) c_longdouble;
pub extern fn __ldexpl(__x: c_longdouble, __exponent: c_int) c_longdouble;
pub extern fn logl(__x: c_longdouble) c_longdouble;
pub extern fn __logl(__x: c_longdouble) c_longdouble;
pub extern fn log10l(__x: c_longdouble) c_longdouble;
pub extern fn __log10l(__x: c_longdouble) c_longdouble;
pub extern fn modfl(__x: c_longdouble, __iptr: [*c]c_longdouble) c_longdouble;
pub extern fn __modfl(__x: c_longdouble, __iptr: [*c]c_longdouble) c_longdouble;
pub extern fn expm1l(__x: c_longdouble) c_longdouble;
pub extern fn __expm1l(__x: c_longdouble) c_longdouble;
pub extern fn log1pl(__x: c_longdouble) c_longdouble;
pub extern fn __log1pl(__x: c_longdouble) c_longdouble;
pub extern fn logbl(__x: c_longdouble) c_longdouble;
pub extern fn __logbl(__x: c_longdouble) c_longdouble;
pub extern fn exp2l(__x: c_longdouble) c_longdouble;
pub extern fn __exp2l(__x: c_longdouble) c_longdouble;
pub extern fn log2l(__x: c_longdouble) c_longdouble;
pub extern fn __log2l(__x: c_longdouble) c_longdouble;
pub extern fn powl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __powl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn sqrtl(__x: c_longdouble) c_longdouble;
pub extern fn __sqrtl(__x: c_longdouble) c_longdouble;
pub extern fn hypotl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __hypotl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn cbrtl(__x: c_longdouble) c_longdouble;
pub extern fn __cbrtl(__x: c_longdouble) c_longdouble;
pub extern fn ceill(__x: c_longdouble) c_longdouble;
pub extern fn fabsl(__x: c_longdouble) c_longdouble;
pub extern fn floorl(__x: c_longdouble) c_longdouble;
pub extern fn fmodl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __fmodl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn isinfl(__value: c_longdouble) c_int;
pub extern fn finitel(__value: c_longdouble) c_int;
pub extern fn dreml(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __dreml(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn significandl(__x: c_longdouble) c_longdouble;
pub extern fn __significandl(__x: c_longdouble) c_longdouble;
pub extern fn copysignl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn nanl(__tagb: [*c]const u8) c_longdouble;
pub extern fn __nanl(__tagb: [*c]const u8) c_longdouble;
pub extern fn isnanl(__value: c_longdouble) c_int;
pub extern fn j0l(c_longdouble) c_longdouble;
pub extern fn __j0l(c_longdouble) c_longdouble;
pub extern fn j1l(c_longdouble) c_longdouble;
pub extern fn __j1l(c_longdouble) c_longdouble;
pub extern fn jnl(c_int, c_longdouble) c_longdouble;
pub extern fn __jnl(c_int, c_longdouble) c_longdouble;
pub extern fn y0l(c_longdouble) c_longdouble;
pub extern fn __y0l(c_longdouble) c_longdouble;
pub extern fn y1l(c_longdouble) c_longdouble;
pub extern fn __y1l(c_longdouble) c_longdouble;
pub extern fn ynl(c_int, c_longdouble) c_longdouble;
pub extern fn __ynl(c_int, c_longdouble) c_longdouble;
pub extern fn erfl(c_longdouble) c_longdouble;
pub extern fn __erfl(c_longdouble) c_longdouble;
pub extern fn erfcl(c_longdouble) c_longdouble;
pub extern fn __erfcl(c_longdouble) c_longdouble;
pub extern fn lgammal(c_longdouble) c_longdouble;
pub extern fn __lgammal(c_longdouble) c_longdouble;
pub extern fn tgammal(c_longdouble) c_longdouble;
pub extern fn __tgammal(c_longdouble) c_longdouble;
pub extern fn gammal(c_longdouble) c_longdouble;
pub extern fn __gammal(c_longdouble) c_longdouble;
pub extern fn lgammal_r(c_longdouble, __signgamp: [*c]c_int) c_longdouble;
pub extern fn __lgammal_r(c_longdouble, __signgamp: [*c]c_int) c_longdouble;
pub extern fn rintl(__x: c_longdouble) c_longdouble;
pub extern fn __rintl(__x: c_longdouble) c_longdouble;
pub extern fn nextafterl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __nextafterl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn nexttowardl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __nexttowardl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn remainderl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __remainderl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn scalbnl(__x: c_longdouble, __n: c_int) c_longdouble;
pub extern fn __scalbnl(__x: c_longdouble, __n: c_int) c_longdouble;
pub extern fn ilogbl(__x: c_longdouble) c_int;
pub extern fn __ilogbl(__x: c_longdouble) c_int;
pub extern fn scalblnl(__x: c_longdouble, __n: c_long) c_longdouble;
pub extern fn __scalblnl(__x: c_longdouble, __n: c_long) c_longdouble;
pub extern fn nearbyintl(__x: c_longdouble) c_longdouble;
pub extern fn __nearbyintl(__x: c_longdouble) c_longdouble;
pub extern fn roundl(__x: c_longdouble) c_longdouble;
pub extern fn truncl(__x: c_longdouble) c_longdouble;
pub extern fn remquol(__x: c_longdouble, __y: c_longdouble, __quo: [*c]c_int) c_longdouble;
pub extern fn __remquol(__x: c_longdouble, __y: c_longdouble, __quo: [*c]c_int) c_longdouble;
pub extern fn lrintl(__x: c_longdouble) c_long;
pub extern fn __lrintl(__x: c_longdouble) c_long;
pub extern fn llrintl(__x: c_longdouble) c_longlong;
pub extern fn __llrintl(__x: c_longdouble) c_longlong;
pub extern fn lroundl(__x: c_longdouble) c_long;
pub extern fn __lroundl(__x: c_longdouble) c_long;
pub extern fn llroundl(__x: c_longdouble) c_longlong;
pub extern fn __llroundl(__x: c_longdouble) c_longlong;
pub extern fn fdiml(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn __fdiml(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn fmaxl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn fminl(__x: c_longdouble, __y: c_longdouble) c_longdouble;
pub extern fn fmal(__x: c_longdouble, __y: c_longdouble, __z: c_longdouble) c_longdouble;
pub extern fn __fmal(__x: c_longdouble, __y: c_longdouble, __z: c_longdouble) c_longdouble;
pub extern fn scalbl(__x: c_longdouble, __n: c_longdouble) c_longdouble;
pub extern fn __scalbl(__x: c_longdouble, __n: c_longdouble) c_longdouble;
pub extern fn __fpclassifyf128(__value: _Float128) c_int;
pub extern fn __signbitf128(__value: _Float128) c_int;
pub extern fn __isinff128(__value: _Float128) c_int;
pub extern fn __finitef128(__value: _Float128) c_int;
pub extern fn __isnanf128(__value: _Float128) c_int;
pub extern fn __iseqsigf128(__x: _Float128, __y: _Float128) c_int;
pub extern fn __issignalingf128(__value: _Float128) c_int;
pub extern var signgam: c_int;
pub const FP_NAN: c_int = 0;
pub const FP_INFINITE: c_int = 1;
pub const FP_ZERO: c_int = 2;
pub const FP_SUBNORMAL: c_int = 3;
pub const FP_NORMAL: c_int = 4;
const enum_unnamed_5 = c_uint;
pub extern fn gsl_log1p(x: f64) f64;
pub extern fn gsl_expm1(x: f64) f64;
pub extern fn gsl_hypot(x: f64, y: f64) f64;
pub extern fn gsl_hypot3(x: f64, y: f64, z: f64) f64;
pub extern fn gsl_acosh(x: f64) f64;
pub extern fn gsl_asinh(x: f64) f64;
pub extern fn gsl_atanh(x: f64) f64;
pub extern fn gsl_isnan(x: f64) c_int;
pub extern fn gsl_isinf(x: f64) c_int;
pub extern fn gsl_finite(x: f64) c_int;
pub extern fn gsl_nan() f64;
pub extern fn gsl_posinf() f64;
pub extern fn gsl_neginf() f64;
pub extern fn gsl_fdiv(x: f64, y: f64) f64;
pub extern fn gsl_coerce_double(x: f64) f64;
pub extern fn gsl_coerce_float(x: f32) f32;
pub extern fn gsl_coerce_long_double(x: c_longdouble) c_longdouble;
pub extern fn gsl_ldexp(x: f64, e: c_int) f64;
pub extern fn gsl_frexp(x: f64, e: [*c]c_int) f64;
pub extern fn gsl_fcmp(x1: f64, x2: f64, epsilon: f64) c_int;
pub const gsl_prec_t = c_uint;
pub const gsl_prec_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_eps",
});
pub const gsl_prec_sqrt_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_sqrt_eps",
});
pub const gsl_prec_root3_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_root3_eps",
});
pub const gsl_prec_root4_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_root4_eps",
});
pub const gsl_prec_root5_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_root5_eps",
});
pub const gsl_prec_root6_eps: [*c]const f64 = @extern([*c]const f64, .{
    .name = "gsl_prec_root6_eps",
});
pub extern fn gsl_pow_2(x: f64) f64;
pub extern fn gsl_pow_3(x: f64) f64;
pub extern fn gsl_pow_4(x: f64) f64;
pub extern fn gsl_pow_5(x: f64) f64;
pub extern fn gsl_pow_6(x: f64) f64;
pub extern fn gsl_pow_7(x: f64) f64;
pub extern fn gsl_pow_8(x: f64) f64;
pub extern fn gsl_pow_9(x: f64) f64;
pub extern fn gsl_pow_int(x: f64, n: c_int) f64;
pub extern fn gsl_pow_uint(x: f64, n: c_uint) f64;
pub extern fn gsl_max(a: f64, b: f64) f64;
pub extern fn gsl_min(a: f64, b: f64) f64;
pub const struct_gsl_function_struct = extern struct {
    function: ?*const fn (f64, ?*anyopaque) callconv(.c) f64 = @import("std").mem.zeroes(?*const fn (f64, ?*anyopaque) callconv(.c) f64),
    params: ?*anyopaque = @import("std").mem.zeroes(?*anyopaque),
};
pub const gsl_function = struct_gsl_function_struct;
pub const struct_gsl_function_fdf_struct = extern struct {
    f: ?*const fn (f64, ?*anyopaque) callconv(.c) f64 = @import("std").mem.zeroes(?*const fn (f64, ?*anyopaque) callconv(.c) f64),
    df: ?*const fn (f64, ?*anyopaque) callconv(.c) f64 = @import("std").mem.zeroes(?*const fn (f64, ?*anyopaque) callconv(.c) f64),
    fdf: ?*const fn (f64, ?*anyopaque, [*c]f64, [*c]f64) callconv(.c) void = @import("std").mem.zeroes(?*const fn (f64, ?*anyopaque, [*c]f64, [*c]f64) callconv(.c) void),
    params: ?*anyopaque = @import("std").mem.zeroes(?*anyopaque),
};
pub const gsl_function_fdf = struct_gsl_function_fdf_struct;
pub const struct_gsl_function_vec_struct = extern struct {
    function: ?*const fn (f64, [*c]f64, ?*anyopaque) callconv(.c) c_int = @import("std").mem.zeroes(?*const fn (f64, [*c]f64, ?*anyopaque) callconv(.c) c_int),
    params: ?*anyopaque = @import("std").mem.zeroes(?*anyopaque),
};
pub const gsl_function_vec = struct_gsl_function_vec_struct;
pub extern fn gsl_blas_sdsdot(alpha: f32, X: [*c]const gsl_vector_float, Y: [*c]const gsl_vector_float, result: [*c]f32) c_int;
pub extern fn gsl_blas_dsdot(X: [*c]const gsl_vector_float, Y: [*c]const gsl_vector_float, result: [*c]f64) c_int;
pub extern fn gsl_blas_sdot(X: [*c]const gsl_vector_float, Y: [*c]const gsl_vector_float, result: [*c]f32) c_int;
pub extern fn gsl_blas_ddot(X: [*c]const gsl_vector, Y: [*c]const gsl_vector, result: [*c]f64) c_int;
pub extern fn gsl_blas_cdotu(X: [*c]const gsl_vector_complex_float, Y: [*c]const gsl_vector_complex_float, dotu: [*c]gsl_complex_float) c_int;
pub extern fn gsl_blas_cdotc(X: [*c]const gsl_vector_complex_float, Y: [*c]const gsl_vector_complex_float, dotc: [*c]gsl_complex_float) c_int;
pub extern fn gsl_blas_zdotu(X: [*c]const gsl_vector_complex, Y: [*c]const gsl_vector_complex, dotu: [*c]gsl_complex) c_int;
pub extern fn gsl_blas_zdotc(X: [*c]const gsl_vector_complex, Y: [*c]const gsl_vector_complex, dotc: [*c]gsl_complex) c_int;
pub extern fn gsl_blas_snrm2(X: [*c]const gsl_vector_float) f32;
pub extern fn gsl_blas_sasum(X: [*c]const gsl_vector_float) f32;
pub extern fn gsl_blas_dnrm2(X: [*c]const gsl_vector) f64;
pub extern fn gsl_blas_dasum(X: [*c]const gsl_vector) f64;
pub extern fn gsl_blas_scnrm2(X: [*c]const gsl_vector_complex_float) f32;
pub extern fn gsl_blas_scasum(X: [*c]const gsl_vector_complex_float) f32;
pub extern fn gsl_blas_dznrm2(X: [*c]const gsl_vector_complex) f64;
pub extern fn gsl_blas_dzasum(X: [*c]const gsl_vector_complex) f64;
pub extern fn gsl_blas_isamax(X: [*c]const gsl_vector_float) CBLAS_INDEX_t;
pub extern fn gsl_blas_idamax(X: [*c]const gsl_vector) CBLAS_INDEX_t;
pub extern fn gsl_blas_icamax(X: [*c]const gsl_vector_complex_float) CBLAS_INDEX_t;
pub extern fn gsl_blas_izamax(X: [*c]const gsl_vector_complex) CBLAS_INDEX_t;
pub extern fn gsl_blas_sswap(X: [*c]gsl_vector_float, Y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_scopy(X: [*c]const gsl_vector_float, Y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_saxpy(alpha: f32, X: [*c]const gsl_vector_float, Y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_dswap(X: [*c]gsl_vector, Y: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_dcopy(X: [*c]const gsl_vector, Y: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_daxpy(alpha: f64, X: [*c]const gsl_vector, Y: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_cswap(X: [*c]gsl_vector_complex_float, Y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_ccopy(X: [*c]const gsl_vector_complex_float, Y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_caxpy(alpha: gsl_complex_float, X: [*c]const gsl_vector_complex_float, Y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_zswap(X: [*c]gsl_vector_complex, Y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_zcopy(X: [*c]const gsl_vector_complex, Y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_zaxpy(alpha: gsl_complex, X: [*c]const gsl_vector_complex, Y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_srotg(a: [*c]f32, b: [*c]f32, c: [*c]f32, s: [*c]f32) c_int;
pub extern fn gsl_blas_srotmg(d1: [*c]f32, d2: [*c]f32, b1: [*c]f32, b2: f32, P: [*c]f32) c_int;
pub extern fn gsl_blas_srot(X: [*c]gsl_vector_float, Y: [*c]gsl_vector_float, c: f32, s: f32) c_int;
pub extern fn gsl_blas_srotm(X: [*c]gsl_vector_float, Y: [*c]gsl_vector_float, P: [*c]const f32) c_int;
pub extern fn gsl_blas_drotg(a: [*c]f64, b: [*c]f64, c: [*c]f64, s: [*c]f64) c_int;
pub extern fn gsl_blas_drotmg(d1: [*c]f64, d2: [*c]f64, b1: [*c]f64, b2: f64, P: [*c]f64) c_int;
pub extern fn gsl_blas_drot(X: [*c]gsl_vector, Y: [*c]gsl_vector, c: f64, s: f64) c_int;
pub extern fn gsl_blas_drotm(X: [*c]gsl_vector, Y: [*c]gsl_vector, P: [*c]const f64) c_int;
pub extern fn gsl_blas_sscal(alpha: f32, X: [*c]gsl_vector_float) void;
pub extern fn gsl_blas_dscal(alpha: f64, X: [*c]gsl_vector) void;
pub extern fn gsl_blas_cscal(alpha: gsl_complex_float, X: [*c]gsl_vector_complex_float) void;
pub extern fn gsl_blas_zscal(alpha: gsl_complex, X: [*c]gsl_vector_complex) void;
pub extern fn gsl_blas_csscal(alpha: f32, X: [*c]gsl_vector_complex_float) void;
pub extern fn gsl_blas_zdscal(alpha: f64, X: [*c]gsl_vector_complex) void;
pub extern fn gsl_blas_sgemv(TransA: CBLAS_TRANSPOSE_t, alpha: f32, A: [*c]const gsl_matrix_float, X: [*c]const gsl_vector_float, beta: f32, Y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_strmv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_float, X: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_strsv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_float, X: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_dgemv(TransA: CBLAS_TRANSPOSE_t, alpha: f64, A: [*c]const gsl_matrix, X: [*c]const gsl_vector, beta: f64, Y: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_dtrmv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix, X: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_dtrsv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix, X: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_cgemv(TransA: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, X: [*c]const gsl_vector_complex_float, beta: gsl_complex_float, Y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_ctrmv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_complex_float, X: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_ctrsv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_complex_float, X: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_zgemv(TransA: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, X: [*c]const gsl_vector_complex, beta: gsl_complex, Y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_ztrmv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_complex, X: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_ztrsv(Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, A: [*c]const gsl_matrix_complex, X: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_ssymv(Uplo: CBLAS_UPLO_t, alpha: f32, A: [*c]const gsl_matrix_float, X: [*c]const gsl_vector_float, beta: f32, Y: [*c]gsl_vector_float) c_int;
pub extern fn gsl_blas_sger(alpha: f32, X: [*c]const gsl_vector_float, Y: [*c]const gsl_vector_float, A: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_ssyr(Uplo: CBLAS_UPLO_t, alpha: f32, X: [*c]const gsl_vector_float, A: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_ssyr2(Uplo: CBLAS_UPLO_t, alpha: f32, X: [*c]const gsl_vector_float, Y: [*c]const gsl_vector_float, A: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_dsymv(Uplo: CBLAS_UPLO_t, alpha: f64, A: [*c]const gsl_matrix, X: [*c]const gsl_vector, beta: f64, Y: [*c]gsl_vector) c_int;
pub extern fn gsl_blas_dger(alpha: f64, X: [*c]const gsl_vector, Y: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dsyr(Uplo: CBLAS_UPLO_t, alpha: f64, X: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dsyr2(Uplo: CBLAS_UPLO_t, alpha: f64, X: [*c]const gsl_vector, Y: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_chemv(Uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, X: [*c]const gsl_vector_complex_float, beta: gsl_complex_float, Y: [*c]gsl_vector_complex_float) c_int;
pub extern fn gsl_blas_cgeru(alpha: gsl_complex_float, X: [*c]const gsl_vector_complex_float, Y: [*c]const gsl_vector_complex_float, A: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_cgerc(alpha: gsl_complex_float, X: [*c]const gsl_vector_complex_float, Y: [*c]const gsl_vector_complex_float, A: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_cher(Uplo: CBLAS_UPLO_t, alpha: f32, X: [*c]const gsl_vector_complex_float, A: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_cher2(Uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, X: [*c]const gsl_vector_complex_float, Y: [*c]const gsl_vector_complex_float, A: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_zhemv(Uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, X: [*c]const gsl_vector_complex, beta: gsl_complex, Y: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_blas_zgeru(alpha: gsl_complex, X: [*c]const gsl_vector_complex, Y: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zgerc(alpha: gsl_complex, X: [*c]const gsl_vector_complex, Y: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zher(Uplo: CBLAS_UPLO_t, alpha: f64, X: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zher2(Uplo: CBLAS_UPLO_t, alpha: gsl_complex, X: [*c]const gsl_vector_complex, Y: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_sgemm(TransA: CBLAS_TRANSPOSE_t, TransB: CBLAS_TRANSPOSE_t, alpha: f32, A: [*c]const gsl_matrix_float, B: [*c]const gsl_matrix_float, beta: f32, C: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_ssymm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: f32, A: [*c]const gsl_matrix_float, B: [*c]const gsl_matrix_float, beta: f32, C: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_ssyrk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f32, A: [*c]const gsl_matrix_float, beta: f32, C: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_ssyr2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f32, A: [*c]const gsl_matrix_float, B: [*c]const gsl_matrix_float, beta: f32, C: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_strmm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: f32, A: [*c]const gsl_matrix_float, B: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_strsm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: f32, A: [*c]const gsl_matrix_float, B: [*c]gsl_matrix_float) c_int;
pub extern fn gsl_blas_dgemm(TransA: CBLAS_TRANSPOSE_t, TransB: CBLAS_TRANSPOSE_t, alpha: f64, A: [*c]const gsl_matrix, B: [*c]const gsl_matrix, beta: f64, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dsymm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: f64, A: [*c]const gsl_matrix, B: [*c]const gsl_matrix, beta: f64, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dsyrk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f64, A: [*c]const gsl_matrix, beta: f64, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dsyr2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f64, A: [*c]const gsl_matrix, B: [*c]const gsl_matrix, beta: f64, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dtrmm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: f64, A: [*c]const gsl_matrix, B: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_dtrsm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: f64, A: [*c]const gsl_matrix, B: [*c]gsl_matrix) c_int;
pub extern fn gsl_blas_cgemm(TransA: CBLAS_TRANSPOSE_t, TransB: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]const gsl_matrix_complex_float, beta: gsl_complex_float, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_csymm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]const gsl_matrix_complex_float, beta: gsl_complex_float, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_csyrk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, beta: gsl_complex_float, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_csyr2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]const gsl_matrix_complex_float, beta: gsl_complex_float, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_ctrmm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_ctrsm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_zgemm(TransA: CBLAS_TRANSPOSE_t, TransB: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, beta: gsl_complex, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zsymm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, beta: gsl_complex, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zsyrk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, beta: gsl_complex, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zsyr2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, beta: gsl_complex, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_ztrmm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_ztrsm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, TransA: CBLAS_TRANSPOSE_t, Diag: CBLAS_DIAG_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_chemm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]const gsl_matrix_complex_float, beta: gsl_complex_float, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_cherk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f32, A: [*c]const gsl_matrix_complex_float, beta: f32, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_cher2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: [*c]const gsl_matrix_complex_float, B: [*c]const gsl_matrix_complex_float, beta: f32, C: [*c]gsl_matrix_complex_float) c_int;
pub extern fn gsl_blas_zhemm(Side: CBLAS_SIDE_t, Uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, beta: gsl_complex, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zherk(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: f64, A: [*c]const gsl_matrix_complex, beta: f64, C: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_blas_zher2k(Uplo: CBLAS_UPLO_t, Trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, beta: f64, C: [*c]gsl_matrix_complex) c_int;
pub const GSL_LINALG_MOD_NONE: c_int = 0;
pub const GSL_LINALG_MOD_TRANSPOSE: c_int = 1;
pub const GSL_LINALG_MOD_CONJUGATE: c_int = 2;
pub const gsl_linalg_matrix_mod_t = c_uint;
pub extern fn gsl_linalg_matmult(A: [*c]const gsl_matrix, B: [*c]const gsl_matrix, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_matmult_mod(A: [*c]const gsl_matrix, modA: gsl_linalg_matrix_mod_t, B: [*c]const gsl_matrix, modB: gsl_linalg_matrix_mod_t, C: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_exponential_ss(A: [*c]const gsl_matrix, eA: [*c]gsl_matrix, mode: gsl_mode_t) c_int;
pub extern fn gsl_linalg_householder_transform(v: [*c]gsl_vector) f64;
pub extern fn gsl_linalg_householder_transform2(alpha: [*c]f64, v: [*c]gsl_vector) f64;
pub extern fn gsl_linalg_complex_householder_transform(v: [*c]gsl_vector_complex) gsl_complex;
pub extern fn gsl_linalg_householder_hm(tau: f64, v: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_householder_mh(tau: f64, v: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_householder_hv(tau: f64, v: [*c]const gsl_vector, w: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_householder_left(tau: f64, v: [*c]const gsl_vector, A: [*c]gsl_matrix, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_householder_right(tau: f64, v: [*c]const gsl_vector, A: [*c]gsl_matrix, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_householder_hm1(tau: f64, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_complex_householder_hm(tau: gsl_complex, v: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_householder_mh(tau: gsl_complex, v: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_householder_hv(tau: gsl_complex, v: [*c]const gsl_vector_complex, w: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_householder_left(tau: gsl_complex, v: [*c]const gsl_vector_complex, A: [*c]gsl_matrix_complex, work: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_hessenberg_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_hessenberg_unpack(H: [*c]gsl_matrix, tau: [*c]gsl_vector, U: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_hessenberg_unpack_accum(H: [*c]gsl_matrix, tau: [*c]gsl_vector, U: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_hessenberg_set_zero(H: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_hessenberg_submatrix(M: [*c]gsl_matrix, A: [*c]gsl_matrix, top: usize, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_hesstri_decomp(A: [*c]gsl_matrix, B: [*c]gsl_matrix, U: [*c]gsl_matrix, V: [*c]gsl_matrix, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_decomp(A: [*c]gsl_matrix, V: [*c]gsl_matrix, S: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_decomp_mod(A: [*c]gsl_matrix, X: [*c]gsl_matrix, V: [*c]gsl_matrix, S: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_decomp_jacobi(A: [*c]gsl_matrix, Q: [*c]gsl_matrix, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_solve(U: [*c]const gsl_matrix, Q: [*c]const gsl_matrix, S: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_solve2(tol: f64, U: [*c]const gsl_matrix, V: [*c]const gsl_matrix, S: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_lssolve(lambda: f64, U: [*c]const gsl_matrix, V: [*c]const gsl_matrix, S: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector, rnorm: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_SV_leverage(U: [*c]const gsl_matrix, h: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_decomp(A: [*c]gsl_matrix, p: [*c]gsl_permutation, signum: [*c]c_int) c_int;
pub extern fn gsl_linalg_LU_solve(LU: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_svx(LU: [*c]const gsl_matrix, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_refine(A: [*c]const gsl_matrix, LU: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_invert(LU: [*c]const gsl_matrix, p: [*c]const gsl_permutation, inverse: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_LU_invx(LU: [*c]gsl_matrix, p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_linalg_LU_det(LU: [*c]gsl_matrix, signum: c_int) f64;
pub extern fn gsl_linalg_LU_lndet(LU: [*c]gsl_matrix) f64;
pub extern fn gsl_linalg_LU_sgndet(lu: [*c]gsl_matrix, signum: c_int) c_int;
pub extern fn gsl_linalg_LU_band_decomp(M: usize, lb: usize, ub: usize, AB: [*c]gsl_matrix, piv: [*c]gsl_vector_uint) c_int;
pub extern fn gsl_linalg_LU_band_solve(lb: usize, ub: usize, LUB: [*c]const gsl_matrix, piv: [*c]const gsl_vector_uint, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_band_svx(lb: usize, ub: usize, LUB: [*c]const gsl_matrix, piv: [*c]const gsl_vector_uint, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LU_band_unpack(M: usize, lb: usize, ub: usize, LUB: [*c]const gsl_matrix, piv: [*c]const gsl_vector_uint, L: [*c]gsl_matrix, U: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_complex_LU_decomp(A: [*c]gsl_matrix_complex, p: [*c]gsl_permutation, signum: [*c]c_int) c_int;
pub extern fn gsl_linalg_complex_LU_solve(LU: [*c]const gsl_matrix_complex, p: [*c]const gsl_permutation, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_LU_svx(LU: [*c]const gsl_matrix_complex, p: [*c]const gsl_permutation, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_LU_refine(A: [*c]const gsl_matrix_complex, LU: [*c]const gsl_matrix_complex, p: [*c]const gsl_permutation, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex, work: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_LU_invert(LU: [*c]const gsl_matrix_complex, p: [*c]const gsl_permutation, inverse: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_LU_invx(LU: [*c]gsl_matrix_complex, p: [*c]const gsl_permutation) c_int;
pub extern fn gsl_linalg_complex_LU_det(LU: [*c]gsl_matrix_complex, signum: c_int) gsl_complex;
pub extern fn gsl_linalg_complex_LU_lndet(LU: [*c]gsl_matrix_complex) f64;
pub extern fn gsl_linalg_complex_LU_sgndet(LU: [*c]gsl_matrix_complex, signum: c_int) gsl_complex;
pub extern fn gsl_linalg_QR_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_decomp_old(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_decomp_r(A: [*c]gsl_matrix, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_solve(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_solve_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_svx(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_lssolve(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_lssolve_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_lssolvem_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, B: [*c]const gsl_matrix, X: [*c]gsl_matrix, work: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_QRsolve(Q: [*c]gsl_matrix, R: [*c]gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_Rsolve(QR: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_Rsvx(QR: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_update(Q: [*c]gsl_matrix, R: [*c]gsl_matrix, w: [*c]gsl_vector, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_linalg_QR_QTvec(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, v: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_QTvec_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_Qvec(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, v: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_QTmat(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_QTmat_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, B: [*c]gsl_matrix, work: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_matQ(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_unpack(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, Q: [*c]gsl_matrix, R: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_unpack_r(QR: [*c]const gsl_matrix, T: [*c]const gsl_matrix, Q: [*c]gsl_matrix, R: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_R_solve(R: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_R_svx(R: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_rcond(QR: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_complex_QR_decomp(A: [*c]gsl_matrix_complex, tau: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_decomp_r(A: [*c]gsl_matrix_complex, T: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_QR_solve(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_solve_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_svx(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_lssolve(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex, residual: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_lssolve_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex, work: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_lssolvem_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, B: [*c]const gsl_matrix_complex, X: [*c]gsl_matrix_complex, work: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_QR_QHvec(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, v: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_QHvec_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, b: [*c]gsl_vector_complex, work: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_QHmat_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, B: [*c]gsl_matrix_complex, work: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_QR_Qvec(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, v: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_QR_unpack(QR: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, Q: [*c]gsl_matrix_complex, R: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_QR_unpack_r(QR: [*c]const gsl_matrix_complex, T: [*c]const gsl_matrix_complex, Q: [*c]gsl_matrix_complex, R: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_QR_band_decomp_L2(M: usize, p: usize, q: usize, AB: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_band_unpack_L2(p: usize, q: usize, QRB: [*c]const gsl_matrix, tau: [*c]const gsl_vector, Q: [*c]gsl_matrix, R: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QRPT_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector, p: [*c]gsl_permutation, signum: [*c]c_int, norm: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_decomp2(A: [*c]const gsl_matrix, q: [*c]gsl_matrix, r: [*c]gsl_matrix, tau: [*c]gsl_vector, p: [*c]gsl_permutation, signum: [*c]c_int, norm: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_solve(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_lssolve(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_lssolve2(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, rank: usize, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_svx(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_QRsolve(Q: [*c]const gsl_matrix, R: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_Rsolve(QR: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_Rsvx(QR: [*c]const gsl_matrix, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_update(Q: [*c]gsl_matrix, R: [*c]gsl_matrix, p: [*c]const gsl_permutation, u: [*c]gsl_vector, v: [*c]const gsl_vector) c_int;
pub extern fn gsl_linalg_QRPT_rank(QR: [*c]const gsl_matrix, tol: f64) usize;
pub extern fn gsl_linalg_QRPT_rcond(QR: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UD_decomp(U: [*c]gsl_matrix, D: [*c]const gsl_vector, Y: [*c]gsl_matrix, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_UD_lssolve(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UD_lssvx(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UD_QTvec(Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UR_decomp(S: [*c]gsl_matrix, A: [*c]gsl_matrix, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_UR_lssolve(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UR_lssvx(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UR_QTvec(Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UU_decomp(U: [*c]gsl_matrix, S: [*c]gsl_matrix, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QR_UU_lssolve(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UU_lssvx(R: [*c]const gsl_matrix, Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, x: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UU_QTvec(Y: [*c]const gsl_matrix, T: [*c]const gsl_matrix, b: [*c]gsl_vector, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QR_UZ_decomp(S: [*c]gsl_matrix, A: [*c]gsl_matrix, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_QL_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_QL_unpack(QL: [*c]const gsl_matrix, tau: [*c]const gsl_vector, Q: [*c]gsl_matrix, L: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_COD_decomp(A: [*c]gsl_matrix, tau_Q: [*c]gsl_vector, tau_Z: [*c]gsl_vector, p: [*c]gsl_permutation, rank: [*c]usize, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_COD_decomp_e(A: [*c]gsl_matrix, tau_Q: [*c]gsl_vector, tau_Z: [*c]gsl_vector, p: [*c]gsl_permutation, tol: f64, rank: [*c]usize, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_COD_lssolve(QRZT: [*c]const gsl_matrix, tau_Q: [*c]const gsl_vector, tau_Z: [*c]const gsl_vector, perm: [*c]const gsl_permutation, rank: usize, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_COD_lssolve2(lambda: f64, QRZT: [*c]const gsl_matrix, tau_Q: [*c]const gsl_vector, tau_Z: [*c]const gsl_vector, perm: [*c]const gsl_permutation, rank: usize, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector, S: [*c]gsl_matrix, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_COD_unpack(QRZT: [*c]const gsl_matrix, tau_Q: [*c]const gsl_vector, tau_Z: [*c]const gsl_vector, rank: usize, Q: [*c]gsl_matrix, R: [*c]gsl_matrix, Z: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_COD_matZ(QRZT: [*c]const gsl_matrix, tau_Z: [*c]const gsl_vector, rank: usize, A: [*c]gsl_matrix, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_lssolve(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_QTvec(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, v: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_solve_T(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_svx_T(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_lssolve_T(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector, residual: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_Lsolve_T(LQ: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_Lsvx_T(LQ: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_L_solve_T(L: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_vecQ(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, v: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_vecQT(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, v: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_unpack(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, Q: [*c]gsl_matrix, L: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_LQ_update(Q: [*c]gsl_matrix, R: [*c]gsl_matrix, v: [*c]const gsl_vector, w: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_LQ_LQsolve(Q: [*c]gsl_matrix, L: [*c]gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector, p: [*c]gsl_permutation, signum: [*c]c_int, norm: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_decomp2(A: [*c]const gsl_matrix, q: [*c]gsl_matrix, r: [*c]gsl_matrix, tau: [*c]gsl_vector, p: [*c]gsl_permutation, signum: [*c]c_int, norm: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_solve_T(QR: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_svx_T(LQ: [*c]const gsl_matrix, tau: [*c]const gsl_vector, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_LQsolve_T(Q: [*c]const gsl_matrix, L: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_Lsolve_T(LQ: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_Lsvx_T(LQ: [*c]const gsl_matrix, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_PTLQ_update(Q: [*c]gsl_matrix, L: [*c]gsl_matrix, p: [*c]const gsl_permutation, v: [*c]const gsl_vector, w: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_decomp(A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_decomp1(A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_solve(cholesky: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_solve_mat(cholesky: [*c]const gsl_matrix, B: [*c]const gsl_matrix, X: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_svx(cholesky: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_svx_mat(cholesky: [*c]const gsl_matrix, X: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_invert(cholesky: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_decomp_unit(A: [*c]gsl_matrix, D: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_scale(A: [*c]const gsl_matrix, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_scale_apply(A: [*c]gsl_matrix, S: [*c]const gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_decomp2(A: [*c]gsl_matrix, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_svx2(LLT: [*c]const gsl_matrix, S: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_solve2(LLT: [*c]const gsl_matrix, S: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_rcond(LLT: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_complex_cholesky_decomp(A: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_cholesky_solve(cholesky: [*c]const gsl_matrix_complex, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_cholesky_svx(cholesky: [*c]const gsl_matrix_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_cholesky_invert(cholesky: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_cholesky_scale(A: [*c]const gsl_matrix_complex, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_complex_cholesky_scale_apply(A: [*c]gsl_matrix_complex, S: [*c]const gsl_vector) c_int;
pub extern fn gsl_linalg_complex_cholesky_decomp2(A: [*c]gsl_matrix_complex, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_complex_cholesky_svx2(LLT: [*c]const gsl_matrix_complex, S: [*c]const gsl_vector, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_complex_cholesky_solve2(LLT: [*c]const gsl_matrix_complex, S: [*c]const gsl_vector, b: [*c]const gsl_vector_complex, x: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_pcholesky_decomp(A: [*c]gsl_matrix, p: [*c]gsl_permutation) c_int;
pub extern fn gsl_linalg_pcholesky_solve(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_pcholesky_svx(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_pcholesky_decomp2(A: [*c]gsl_matrix, p: [*c]gsl_permutation, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_pcholesky_solve2(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, S: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_pcholesky_svx2(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, S: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_pcholesky_invert(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, Ainv: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_pcholesky_rcond(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_mcholesky_decomp(A: [*c]gsl_matrix, p: [*c]gsl_permutation, E: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_mcholesky_solve(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_mcholesky_svx(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_mcholesky_rcond(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_mcholesky_invert(LDLT: [*c]const gsl_matrix, p: [*c]const gsl_permutation, Ainv: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_decomp(A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_solve(LLT: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_band_svx(LLT: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_band_solvem(LLT: [*c]const gsl_matrix, B: [*c]const gsl_matrix, X: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_svxm(LLT: [*c]const gsl_matrix, X: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_invert(LLT: [*c]const gsl_matrix, Ainv: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_unpack(LLT: [*c]const gsl_matrix, L: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_cholesky_band_scale(A: [*c]const gsl_matrix, S: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_band_scale_apply(A: [*c]gsl_matrix, S: [*c]const gsl_vector) c_int;
pub extern fn gsl_linalg_cholesky_band_rcond(LLT: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_decomp(A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_ldlt_solve(LDLT: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_svx(LDLT: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_rcond(LDLT: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_band_decomp(A: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_ldlt_band_solve(LDLT: [*c]const gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_band_svx(LDLT: [*c]const gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_band_unpack(LDLT: [*c]const gsl_matrix, L: [*c]gsl_matrix, D: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_ldlt_band_rcond(LDLT: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_symmtd_decomp(A: [*c]gsl_matrix, tau: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_symmtd_unpack(A: [*c]const gsl_matrix, tau: [*c]const gsl_vector, Q: [*c]gsl_matrix, diag: [*c]gsl_vector, subdiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_symmtd_unpack_T(A: [*c]const gsl_matrix, diag: [*c]gsl_vector, subdiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_hermtd_decomp(A: [*c]gsl_matrix_complex, tau: [*c]gsl_vector_complex) c_int;
pub extern fn gsl_linalg_hermtd_unpack(A: [*c]const gsl_matrix_complex, tau: [*c]const gsl_vector_complex, U: [*c]gsl_matrix_complex, diag: [*c]gsl_vector, sudiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_hermtd_unpack_T(A: [*c]const gsl_matrix_complex, diag: [*c]gsl_vector, subdiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_HH_solve(A: [*c]gsl_matrix, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_HH_svx(A: [*c]gsl_matrix, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_solve_symm_tridiag(diag: [*c]const gsl_vector, offdiag: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_solve_tridiag(diag: [*c]const gsl_vector, abovediag: [*c]const gsl_vector, belowdiag: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_solve_symm_cyc_tridiag(diag: [*c]const gsl_vector, offdiag: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_solve_cyc_tridiag(diag: [*c]const gsl_vector, abovediag: [*c]const gsl_vector, belowdiag: [*c]const gsl_vector, b: [*c]const gsl_vector, x: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_bidiag_decomp(A: [*c]gsl_matrix, tau_U: [*c]gsl_vector, tau_V: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_bidiag_unpack(A: [*c]const gsl_matrix, tau_U: [*c]const gsl_vector, U: [*c]gsl_matrix, tau_V: [*c]const gsl_vector, V: [*c]gsl_matrix, diag: [*c]gsl_vector, superdiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_bidiag_unpack2(A: [*c]gsl_matrix, tau_U: [*c]gsl_vector, tau_V: [*c]gsl_vector, V: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_bidiag_unpack_B(A: [*c]const gsl_matrix, diag: [*c]gsl_vector, superdiag: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_balance_matrix(A: [*c]gsl_matrix, D: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_balance_accum(A: [*c]gsl_matrix, D: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_balance_columns(A: [*c]gsl_matrix, D: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_tri_rcond(Uplo: CBLAS_UPLO_t, A: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_tri_upper_rcond(A: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_tri_lower_rcond(A: [*c]const gsl_matrix, rcond: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_invnorm1(N: usize, Ainvx: ?*const fn (CBLAS_TRANSPOSE_t, [*c]gsl_vector, ?*anyopaque) callconv(.c) c_int, params: ?*anyopaque, Ainvnorm: [*c]f64, work: [*c]gsl_vector) c_int;
pub extern fn gsl_linalg_tri_upper_invert(T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_tri_lower_invert(T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_tri_upper_unit_invert(T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_tri_lower_unit_invert(T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_tri_invert(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, T: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_complex_tri_invert(Uplo: CBLAS_UPLO_t, Diag: CBLAS_DIAG_t, T: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_tri_LTL(L: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_tri_UL(LU: [*c]gsl_matrix) c_int;
pub extern fn gsl_linalg_complex_tri_LHL(L: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_complex_tri_UL(LU: [*c]gsl_matrix_complex) c_int;
pub extern fn gsl_linalg_givens(a: f64, b: f64, c: [*c]f64, s: [*c]f64) void;
pub extern fn gsl_linalg_givens_gv(v: [*c]gsl_vector, i: usize, j: usize, c: f64, s: f64) void;
pub extern fn gsl_complex_polar(r: f64, theta: f64) gsl_complex;
pub extern fn gsl_complex_rect(x: f64, y: f64) gsl_complex;
pub extern fn gsl_complex_arg(z: gsl_complex) f64;
pub extern fn gsl_complex_abs(z: gsl_complex) f64;
pub extern fn gsl_complex_abs2(z: gsl_complex) f64;
pub extern fn gsl_complex_logabs(z: gsl_complex) f64;
pub extern fn gsl_complex_add(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sub(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_mul(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_div(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_add_real(a: gsl_complex, x: f64) gsl_complex;
pub extern fn gsl_complex_sub_real(a: gsl_complex, x: f64) gsl_complex;
pub extern fn gsl_complex_mul_real(a: gsl_complex, x: f64) gsl_complex;
pub extern fn gsl_complex_div_real(a: gsl_complex, x: f64) gsl_complex;
pub extern fn gsl_complex_add_imag(a: gsl_complex, y: f64) gsl_complex;
pub extern fn gsl_complex_sub_imag(a: gsl_complex, y: f64) gsl_complex;
pub extern fn gsl_complex_mul_imag(a: gsl_complex, y: f64) gsl_complex;
pub extern fn gsl_complex_div_imag(a: gsl_complex, y: f64) gsl_complex;
pub extern fn gsl_complex_conjugate(z: gsl_complex) gsl_complex;
pub extern fn gsl_complex_inverse(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_negative(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sqrt(z: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sqrt_real(x: f64) gsl_complex;
pub extern fn gsl_complex_pow(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_pow_real(a: gsl_complex, b: f64) gsl_complex;
pub extern fn gsl_complex_exp(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_log(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_log10(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_log_b(a: gsl_complex, b: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sin(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_cos(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sec(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_csc(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_tan(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_cot(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arcsin(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arcsin_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arccos(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccos_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arcsec(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arcsec_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arccsc(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccsc_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arctan(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccot(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sinh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_cosh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_sech(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_csch(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_tanh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_coth(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arcsinh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccosh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccosh_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arcsech(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arccsch(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arctanh(a: gsl_complex) gsl_complex;
pub extern fn gsl_complex_arctanh_real(a: f64) gsl_complex;
pub extern fn gsl_complex_arccoth(a: gsl_complex) gsl_complex;
pub const gsl_integration_workspace = extern struct {
    limit: usize = @import("std").mem.zeroes(usize),
    size: usize = @import("std").mem.zeroes(usize),
    nrmax: usize = @import("std").mem.zeroes(usize),
    i: usize = @import("std").mem.zeroes(usize),
    maximum_level: usize = @import("std").mem.zeroes(usize),
    alist: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    blist: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    rlist: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    elist: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    order: [*c]usize = @import("std").mem.zeroes([*c]usize),
    level: [*c]usize = @import("std").mem.zeroes([*c]usize),
};
pub extern fn gsl_integration_workspace_alloc(n: usize) [*c]gsl_integration_workspace;
pub extern fn gsl_integration_workspace_free(w: [*c]gsl_integration_workspace) void;
pub const gsl_integration_qaws_table = extern struct {
    alpha: f64 = @import("std").mem.zeroes(f64),
    beta: f64 = @import("std").mem.zeroes(f64),
    mu: c_int = @import("std").mem.zeroes(c_int),
    nu: c_int = @import("std").mem.zeroes(c_int),
    ri: [25]f64 = @import("std").mem.zeroes([25]f64),
    rj: [25]f64 = @import("std").mem.zeroes([25]f64),
    rg: [25]f64 = @import("std").mem.zeroes([25]f64),
    rh: [25]f64 = @import("std").mem.zeroes([25]f64),
};
pub extern fn gsl_integration_qaws_table_alloc(alpha: f64, beta: f64, mu: c_int, nu: c_int) [*c]gsl_integration_qaws_table;
pub extern fn gsl_integration_qaws_table_set(t: [*c]gsl_integration_qaws_table, alpha: f64, beta: f64, mu: c_int, nu: c_int) c_int;
pub extern fn gsl_integration_qaws_table_free(t: [*c]gsl_integration_qaws_table) void;
pub const GSL_INTEG_COSINE: c_int = 0;
pub const GSL_INTEG_SINE: c_int = 1;
pub const enum_gsl_integration_qawo_enum = c_uint;
pub const gsl_integration_qawo_table = extern struct {
    n: usize = @import("std").mem.zeroes(usize),
    omega: f64 = @import("std").mem.zeroes(f64),
    L: f64 = @import("std").mem.zeroes(f64),
    par: f64 = @import("std").mem.zeroes(f64),
    sine: enum_gsl_integration_qawo_enum = @import("std").mem.zeroes(enum_gsl_integration_qawo_enum),
    chebmo: [*c]f64 = @import("std").mem.zeroes([*c]f64),
};
pub extern fn gsl_integration_qawo_table_alloc(omega: f64, L: f64, sine: enum_gsl_integration_qawo_enum, n: usize) [*c]gsl_integration_qawo_table;
pub extern fn gsl_integration_qawo_table_set(t: [*c]gsl_integration_qawo_table, omega: f64, L: f64, sine: enum_gsl_integration_qawo_enum) c_int;
pub extern fn gsl_integration_qawo_table_set_length(t: [*c]gsl_integration_qawo_table, L: f64) c_int;
pub extern fn gsl_integration_qawo_table_free(t: [*c]gsl_integration_qawo_table) void;
pub const gsl_integration_rule = fn ([*c]const gsl_function, f64, f64, [*c]f64, [*c]f64, [*c]f64, [*c]f64) callconv(.c) void;
pub extern fn gsl_integration_qk15(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qk21(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qk31(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qk41(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qk51(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qk61(f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qcheb(f: [*c]gsl_function, a: f64, b: f64, cheb12: [*c]f64, cheb24: [*c]f64) void;
pub const GSL_INTEG_GAUSS15: c_int = 1;
pub const GSL_INTEG_GAUSS21: c_int = 2;
pub const GSL_INTEG_GAUSS31: c_int = 3;
pub const GSL_INTEG_GAUSS41: c_int = 4;
pub const GSL_INTEG_GAUSS51: c_int = 5;
pub const GSL_INTEG_GAUSS61: c_int = 6;
const enum_unnamed_6 = c_uint;
pub extern fn gsl_integration_qk(n: c_int, xgk: [*c]const f64, wg: [*c]const f64, wgk: [*c]const f64, fv1: [*c]f64, fv2: [*c]f64, f: [*c]const gsl_function, a: f64, b: f64, result: [*c]f64, abserr: [*c]f64, resabs: [*c]f64, resasc: [*c]f64) void;
pub extern fn gsl_integration_qng(f: [*c]const gsl_function, a: f64, b: f64, epsabs: f64, epsrel: f64, result: [*c]f64, abserr: [*c]f64, neval: [*c]usize) c_int;
pub extern fn gsl_integration_qag(f: [*c]const gsl_function, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: usize, key: c_int, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qagi(f: [*c]gsl_function, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qagiu(f: [*c]gsl_function, a: f64, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qagil(f: [*c]gsl_function, b: f64, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qags(f: [*c]const gsl_function, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qagp(f: [*c]const gsl_function, pts: [*c]f64, npts: usize, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qawc(f: [*c]gsl_function, a: f64, b: f64, c: f64, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qaws(f: [*c]gsl_function, a: f64, b: f64, t: [*c]gsl_integration_qaws_table, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qawo(f: [*c]gsl_function, a: f64, epsabs: f64, epsrel: f64, limit: usize, workspace: [*c]gsl_integration_workspace, wf: [*c]gsl_integration_qawo_table, result: [*c]f64, abserr: [*c]f64) c_int;
pub extern fn gsl_integration_qawf(f: [*c]gsl_function, a: f64, epsabs: f64, limit: usize, workspace: [*c]gsl_integration_workspace, cycle_workspace: [*c]gsl_integration_workspace, wf: [*c]gsl_integration_qawo_table, result: [*c]f64, abserr: [*c]f64) c_int;
pub const gsl_integration_glfixed_table = extern struct {
    n: usize = @import("std").mem.zeroes(usize),
    x: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    w: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    precomputed: c_int = @import("std").mem.zeroes(c_int),
};
pub extern fn gsl_integration_glfixed_table_alloc(n: usize) [*c]gsl_integration_glfixed_table;
pub extern fn gsl_integration_glfixed_table_free(t: [*c]gsl_integration_glfixed_table) void;
pub extern fn gsl_integration_glfixed(f: [*c]const gsl_function, a: f64, b: f64, t: [*c]const gsl_integration_glfixed_table) f64;
pub extern fn gsl_integration_glfixed_point(a: f64, b: f64, i: usize, xi: [*c]f64, wi: [*c]f64, t: [*c]const gsl_integration_glfixed_table) c_int;
pub const gsl_integration_cquad_ival = extern struct {
    a: f64 = @import("std").mem.zeroes(f64),
    b: f64 = @import("std").mem.zeroes(f64),
    c: [64]f64 = @import("std").mem.zeroes([64]f64),
    fx: [33]f64 = @import("std").mem.zeroes([33]f64),
    igral: f64 = @import("std").mem.zeroes(f64),
    err: f64 = @import("std").mem.zeroes(f64),
    depth: c_int = @import("std").mem.zeroes(c_int),
    rdepth: c_int = @import("std").mem.zeroes(c_int),
    ndiv: c_int = @import("std").mem.zeroes(c_int),
};
pub const gsl_integration_cquad_workspace = extern struct {
    size: usize = @import("std").mem.zeroes(usize),
    ivals: [*c]gsl_integration_cquad_ival = @import("std").mem.zeroes([*c]gsl_integration_cquad_ival),
    heap: [*c]usize = @import("std").mem.zeroes([*c]usize),
};
pub extern fn gsl_integration_cquad_workspace_alloc(n: usize) [*c]gsl_integration_cquad_workspace;
pub extern fn gsl_integration_cquad_workspace_free(w: [*c]gsl_integration_cquad_workspace) void;
pub extern fn gsl_integration_cquad(f: [*c]const gsl_function, a: f64, b: f64, epsabs: f64, epsrel: f64, ws: [*c]gsl_integration_cquad_workspace, result: [*c]f64, abserr: [*c]f64, nevals: [*c]usize) c_int;
pub const gsl_integration_romberg_workspace = extern struct {
    n: usize = @import("std").mem.zeroes(usize),
    work1: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    work2: [*c]f64 = @import("std").mem.zeroes([*c]f64),
};
pub extern fn gsl_integration_romberg_alloc(n: usize) [*c]gsl_integration_romberg_workspace;
pub extern fn gsl_integration_romberg_free(w: [*c]gsl_integration_romberg_workspace) void;
pub extern fn gsl_integration_romberg(f: [*c]const gsl_function, a: f64, b: f64, epsabs: f64, epsrel: f64, result: [*c]f64, neval: [*c]usize, w: [*c]gsl_integration_romberg_workspace) c_int;
pub const gsl_integration_fixed_params = extern struct {
    alpha: f64 = @import("std").mem.zeroes(f64),
    beta: f64 = @import("std").mem.zeroes(f64),
    a: f64 = @import("std").mem.zeroes(f64),
    b: f64 = @import("std").mem.zeroes(f64),
    zemu: f64 = @import("std").mem.zeroes(f64),
    shft: f64 = @import("std").mem.zeroes(f64),
    slp: f64 = @import("std").mem.zeroes(f64),
    al: f64 = @import("std").mem.zeroes(f64),
    be: f64 = @import("std").mem.zeroes(f64),
};
pub const gsl_integration_fixed_type = extern struct {
    check: ?*const fn (usize, [*c]const gsl_integration_fixed_params) callconv(.c) c_int = @import("std").mem.zeroes(?*const fn (usize, [*c]const gsl_integration_fixed_params) callconv(.c) c_int),
    init: ?*const fn (usize, [*c]f64, [*c]f64, [*c]gsl_integration_fixed_params) callconv(.c) c_int = @import("std").mem.zeroes(?*const fn (usize, [*c]f64, [*c]f64, [*c]gsl_integration_fixed_params) callconv(.c) c_int),
};
pub const gsl_integration_fixed_workspace = extern struct {
    n: usize = @import("std").mem.zeroes(usize),
    weights: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    x: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    diag: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    subdiag: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    type: [*c]const gsl_integration_fixed_type = @import("std").mem.zeroes([*c]const gsl_integration_fixed_type),
};
pub extern var gsl_integration_fixed_legendre: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_chebyshev: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_gegenbauer: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_jacobi: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_laguerre: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_hermite: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_exponential: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_rational: [*c]const gsl_integration_fixed_type;
pub extern var gsl_integration_fixed_chebyshev2: [*c]const gsl_integration_fixed_type;
pub extern fn gsl_integration_fixed_alloc(@"type": [*c]const gsl_integration_fixed_type, n: usize, a: f64, b: f64, alpha: f64, beta: f64) [*c]gsl_integration_fixed_workspace;
pub extern fn gsl_integration_fixed_free(w: [*c]gsl_integration_fixed_workspace) void;
pub extern fn gsl_integration_fixed_n(w: [*c]const gsl_integration_fixed_workspace) usize;
pub extern fn gsl_integration_fixed_nodes(w: [*c]const gsl_integration_fixed_workspace) [*c]f64;
pub extern fn gsl_integration_fixed_weights(w: [*c]const gsl_integration_fixed_workspace) [*c]f64;
pub extern fn gsl_integration_fixed(func: [*c]const gsl_function, result: [*c]f64, w: [*c]const gsl_integration_fixed_workspace) c_int;
pub const gsl_integration_lebedev_workspace = extern struct {
    n: usize = @import("std").mem.zeroes(usize),
    weights: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    x: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    y: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    z: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    theta: [*c]f64 = @import("std").mem.zeroes([*c]f64),
    phi: [*c]f64 = @import("std").mem.zeroes([*c]f64),
};
pub extern fn gsl_integration_lebedev_alloc(n: usize) [*c]gsl_integration_lebedev_workspace;
pub extern fn gsl_integration_lebedev_free(w: [*c]gsl_integration_lebedev_workspace) void;
pub extern fn gsl_integration_lebedev_n(w: [*c]const gsl_integration_lebedev_workspace) usize;
pub const __llvm__ = @as(c_int, 1);
pub const __clang__ = @as(c_int, 1);
pub const __clang_major__ = @as(c_int, 19);
pub const __clang_minor__ = @as(c_int, 1);
pub const __clang_patchlevel__ = @as(c_int, 7);
pub const __clang_version__ = "19.1.7 (https://github.com/ziglang/zig-bootstrap 1c3c59435891bc9caf8cd1d3783773369d191c5f)";
pub const __GNUC__ = @as(c_int, 4);
pub const __GNUC_MINOR__ = @as(c_int, 2);
pub const __GNUC_PATCHLEVEL__ = @as(c_int, 1);
pub const __GXX_ABI_VERSION = @as(c_int, 1002);
pub const __ATOMIC_RELAXED = @as(c_int, 0);
pub const __ATOMIC_CONSUME = @as(c_int, 1);
pub const __ATOMIC_ACQUIRE = @as(c_int, 2);
pub const __ATOMIC_RELEASE = @as(c_int, 3);
pub const __ATOMIC_ACQ_REL = @as(c_int, 4);
pub const __ATOMIC_SEQ_CST = @as(c_int, 5);
pub const __MEMORY_SCOPE_SYSTEM = @as(c_int, 0);
pub const __MEMORY_SCOPE_DEVICE = @as(c_int, 1);
pub const __MEMORY_SCOPE_WRKGRP = @as(c_int, 2);
pub const __MEMORY_SCOPE_WVFRNT = @as(c_int, 3);
pub const __MEMORY_SCOPE_SINGLE = @as(c_int, 4);
pub const __OPENCL_MEMORY_SCOPE_WORK_ITEM = @as(c_int, 0);
pub const __OPENCL_MEMORY_SCOPE_WORK_GROUP = @as(c_int, 1);
pub const __OPENCL_MEMORY_SCOPE_DEVICE = @as(c_int, 2);
pub const __OPENCL_MEMORY_SCOPE_ALL_SVM_DEVICES = @as(c_int, 3);
pub const __OPENCL_MEMORY_SCOPE_SUB_GROUP = @as(c_int, 4);
pub const __FPCLASS_SNAN = @as(c_int, 0x0001);
pub const __FPCLASS_QNAN = @as(c_int, 0x0002);
pub const __FPCLASS_NEGINF = @as(c_int, 0x0004);
pub const __FPCLASS_NEGNORMAL = @as(c_int, 0x0008);
pub const __FPCLASS_NEGSUBNORMAL = @as(c_int, 0x0010);
pub const __FPCLASS_NEGZERO = @as(c_int, 0x0020);
pub const __FPCLASS_POSZERO = @as(c_int, 0x0040);
pub const __FPCLASS_POSSUBNORMAL = @as(c_int, 0x0080);
pub const __FPCLASS_POSNORMAL = @as(c_int, 0x0100);
pub const __FPCLASS_POSINF = @as(c_int, 0x0200);
pub const __PRAGMA_REDEFINE_EXTNAME = @as(c_int, 1);
pub const __VERSION__ = "Clang 19.1.7 (https://github.com/ziglang/zig-bootstrap 1c3c59435891bc9caf8cd1d3783773369d191c5f)";
pub const __OBJC_BOOL_IS_BOOL = @as(c_int, 0);
pub const __CONSTANT_CFSTRINGS__ = @as(c_int, 1);
pub const __clang_literal_encoding__ = "UTF-8";
pub const __clang_wide_literal_encoding__ = "UTF-32";
pub const __ORDER_LITTLE_ENDIAN__ = @as(c_int, 1234);
pub const __ORDER_BIG_ENDIAN__ = @as(c_int, 4321);
pub const __ORDER_PDP_ENDIAN__ = @as(c_int, 3412);
pub const __BYTE_ORDER__ = __ORDER_LITTLE_ENDIAN__;
pub const __LITTLE_ENDIAN__ = @as(c_int, 1);
pub const _LP64 = @as(c_int, 1);
pub const __LP64__ = @as(c_int, 1);
pub const __CHAR_BIT__ = @as(c_int, 8);
pub const __BOOL_WIDTH__ = @as(c_int, 8);
pub const __SHRT_WIDTH__ = @as(c_int, 16);
pub const __INT_WIDTH__ = @as(c_int, 32);
pub const __LONG_WIDTH__ = @as(c_int, 64);
pub const __LLONG_WIDTH__ = @as(c_int, 64);
pub const __BITINT_MAXWIDTH__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 8388608, .decimal);
pub const __SCHAR_MAX__ = @as(c_int, 127);
pub const __SHRT_MAX__ = @as(c_int, 32767);
pub const __INT_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __LONG_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __LONG_LONG_MAX__ = @as(c_longlong, 9223372036854775807);
pub const __WCHAR_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __WCHAR_WIDTH__ = @as(c_int, 32);
pub const __WINT_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_uint, 4294967295, .decimal);
pub const __WINT_WIDTH__ = @as(c_int, 32);
pub const __INTMAX_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __INTMAX_WIDTH__ = @as(c_int, 64);
pub const __SIZE_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __SIZE_WIDTH__ = @as(c_int, 64);
pub const __UINTMAX_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __UINTMAX_WIDTH__ = @as(c_int, 64);
pub const __PTRDIFF_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __PTRDIFF_WIDTH__ = @as(c_int, 64);
pub const __INTPTR_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __INTPTR_WIDTH__ = @as(c_int, 64);
pub const __UINTPTR_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __UINTPTR_WIDTH__ = @as(c_int, 64);
pub const __SIZEOF_DOUBLE__ = @as(c_int, 8);
pub const __SIZEOF_FLOAT__ = @as(c_int, 4);
pub const __SIZEOF_INT__ = @as(c_int, 4);
pub const __SIZEOF_LONG__ = @as(c_int, 8);
pub const __SIZEOF_LONG_DOUBLE__ = @as(c_int, 16);
pub const __SIZEOF_LONG_LONG__ = @as(c_int, 8);
pub const __SIZEOF_POINTER__ = @as(c_int, 8);
pub const __SIZEOF_SHORT__ = @as(c_int, 2);
pub const __SIZEOF_PTRDIFF_T__ = @as(c_int, 8);
pub const __SIZEOF_SIZE_T__ = @as(c_int, 8);
pub const __SIZEOF_WCHAR_T__ = @as(c_int, 4);
pub const __SIZEOF_WINT_T__ = @as(c_int, 4);
pub const __SIZEOF_INT128__ = @as(c_int, 16);
pub const __INTMAX_TYPE__ = c_long;
pub const __INTMAX_FMTd__ = "ld";
pub const __INTMAX_FMTi__ = "li";
pub const __INTMAX_C_SUFFIX__ = @compileError("unable to translate macro: undefined identifier `L`");
// (no file):95:9
pub const __UINTMAX_TYPE__ = c_ulong;
pub const __UINTMAX_FMTo__ = "lo";
pub const __UINTMAX_FMTu__ = "lu";
pub const __UINTMAX_FMTx__ = "lx";
pub const __UINTMAX_FMTX__ = "lX";
pub const __UINTMAX_C_SUFFIX__ = @compileError("unable to translate macro: undefined identifier `UL`");
// (no file):101:9
pub const __PTRDIFF_TYPE__ = c_long;
pub const __PTRDIFF_FMTd__ = "ld";
pub const __PTRDIFF_FMTi__ = "li";
pub const __INTPTR_TYPE__ = c_long;
pub const __INTPTR_FMTd__ = "ld";
pub const __INTPTR_FMTi__ = "li";
pub const __SIZE_TYPE__ = c_ulong;
pub const __SIZE_FMTo__ = "lo";
pub const __SIZE_FMTu__ = "lu";
pub const __SIZE_FMTx__ = "lx";
pub const __SIZE_FMTX__ = "lX";
pub const __WCHAR_TYPE__ = c_int;
pub const __WINT_TYPE__ = c_uint;
pub const __SIG_ATOMIC_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __SIG_ATOMIC_WIDTH__ = @as(c_int, 32);
pub const __CHAR16_TYPE__ = c_ushort;
pub const __CHAR32_TYPE__ = c_uint;
pub const __UINTPTR_TYPE__ = c_ulong;
pub const __UINTPTR_FMTo__ = "lo";
pub const __UINTPTR_FMTu__ = "lu";
pub const __UINTPTR_FMTx__ = "lx";
pub const __UINTPTR_FMTX__ = "lX";
pub const __FLT16_DENORM_MIN__ = @as(f16, 5.9604644775390625e-8);
pub const __FLT16_NORM_MAX__ = @as(f16, 6.5504e+4);
pub const __FLT16_HAS_DENORM__ = @as(c_int, 1);
pub const __FLT16_DIG__ = @as(c_int, 3);
pub const __FLT16_DECIMAL_DIG__ = @as(c_int, 5);
pub const __FLT16_EPSILON__ = @as(f16, 9.765625e-4);
pub const __FLT16_HAS_INFINITY__ = @as(c_int, 1);
pub const __FLT16_HAS_QUIET_NAN__ = @as(c_int, 1);
pub const __FLT16_MANT_DIG__ = @as(c_int, 11);
pub const __FLT16_MAX_10_EXP__ = @as(c_int, 4);
pub const __FLT16_MAX_EXP__ = @as(c_int, 16);
pub const __FLT16_MAX__ = @as(f16, 6.5504e+4);
pub const __FLT16_MIN_10_EXP__ = -@as(c_int, 4);
pub const __FLT16_MIN_EXP__ = -@as(c_int, 13);
pub const __FLT16_MIN__ = @as(f16, 6.103515625e-5);
pub const __FLT_DENORM_MIN__ = @as(f32, 1.40129846e-45);
pub const __FLT_NORM_MAX__ = @as(f32, 3.40282347e+38);
pub const __FLT_HAS_DENORM__ = @as(c_int, 1);
pub const __FLT_DIG__ = @as(c_int, 6);
pub const __FLT_DECIMAL_DIG__ = @as(c_int, 9);
pub const __FLT_EPSILON__ = @as(f32, 1.19209290e-7);
pub const __FLT_HAS_INFINITY__ = @as(c_int, 1);
pub const __FLT_HAS_QUIET_NAN__ = @as(c_int, 1);
pub const __FLT_MANT_DIG__ = @as(c_int, 24);
pub const __FLT_MAX_10_EXP__ = @as(c_int, 38);
pub const __FLT_MAX_EXP__ = @as(c_int, 128);
pub const __FLT_MAX__ = @as(f32, 3.40282347e+38);
pub const __FLT_MIN_10_EXP__ = -@as(c_int, 37);
pub const __FLT_MIN_EXP__ = -@as(c_int, 125);
pub const __FLT_MIN__ = @as(f32, 1.17549435e-38);
pub const __DBL_DENORM_MIN__ = @as(f64, 4.9406564584124654e-324);
pub const __DBL_NORM_MAX__ = @as(f64, 1.7976931348623157e+308);
pub const __DBL_HAS_DENORM__ = @as(c_int, 1);
pub const __DBL_DIG__ = @as(c_int, 15);
pub const __DBL_DECIMAL_DIG__ = @as(c_int, 17);
pub const __DBL_EPSILON__ = @as(f64, 2.2204460492503131e-16);
pub const __DBL_HAS_INFINITY__ = @as(c_int, 1);
pub const __DBL_HAS_QUIET_NAN__ = @as(c_int, 1);
pub const __DBL_MANT_DIG__ = @as(c_int, 53);
pub const __DBL_MAX_10_EXP__ = @as(c_int, 308);
pub const __DBL_MAX_EXP__ = @as(c_int, 1024);
pub const __DBL_MAX__ = @as(f64, 1.7976931348623157e+308);
pub const __DBL_MIN_10_EXP__ = -@as(c_int, 307);
pub const __DBL_MIN_EXP__ = -@as(c_int, 1021);
pub const __DBL_MIN__ = @as(f64, 2.2250738585072014e-308);
pub const __LDBL_DENORM_MIN__ = @as(c_longdouble, 3.64519953188247460253e-4951);
pub const __LDBL_NORM_MAX__ = @as(c_longdouble, 1.18973149535723176502e+4932);
pub const __LDBL_HAS_DENORM__ = @as(c_int, 1);
pub const __LDBL_DIG__ = @as(c_int, 18);
pub const __LDBL_DECIMAL_DIG__ = @as(c_int, 21);
pub const __LDBL_EPSILON__ = @as(c_longdouble, 1.08420217248550443401e-19);
pub const __LDBL_HAS_INFINITY__ = @as(c_int, 1);
pub const __LDBL_HAS_QUIET_NAN__ = @as(c_int, 1);
pub const __LDBL_MANT_DIG__ = @as(c_int, 64);
pub const __LDBL_MAX_10_EXP__ = @as(c_int, 4932);
pub const __LDBL_MAX_EXP__ = @as(c_int, 16384);
pub const __LDBL_MAX__ = @as(c_longdouble, 1.18973149535723176502e+4932);
pub const __LDBL_MIN_10_EXP__ = -@as(c_int, 4931);
pub const __LDBL_MIN_EXP__ = -@as(c_int, 16381);
pub const __LDBL_MIN__ = @as(c_longdouble, 3.36210314311209350626e-4932);
pub const __POINTER_WIDTH__ = @as(c_int, 64);
pub const __BIGGEST_ALIGNMENT__ = @as(c_int, 16);
pub const __WINT_UNSIGNED__ = @as(c_int, 1);
pub const __INT8_TYPE__ = i8;
pub const __INT8_FMTd__ = "hhd";
pub const __INT8_FMTi__ = "hhi";
pub const __INT8_C_SUFFIX__ = "";
pub const __INT16_TYPE__ = c_short;
pub const __INT16_FMTd__ = "hd";
pub const __INT16_FMTi__ = "hi";
pub const __INT16_C_SUFFIX__ = "";
pub const __INT32_TYPE__ = c_int;
pub const __INT32_FMTd__ = "d";
pub const __INT32_FMTi__ = "i";
pub const __INT32_C_SUFFIX__ = "";
pub const __INT64_TYPE__ = c_long;
pub const __INT64_FMTd__ = "ld";
pub const __INT64_FMTi__ = "li";
pub const __INT64_C_SUFFIX__ = @compileError("unable to translate macro: undefined identifier `L`");
// (no file):202:9
pub const __UINT8_TYPE__ = u8;
pub const __UINT8_FMTo__ = "hho";
pub const __UINT8_FMTu__ = "hhu";
pub const __UINT8_FMTx__ = "hhx";
pub const __UINT8_FMTX__ = "hhX";
pub const __UINT8_C_SUFFIX__ = "";
pub const __UINT8_MAX__ = @as(c_int, 255);
pub const __INT8_MAX__ = @as(c_int, 127);
pub const __UINT16_TYPE__ = c_ushort;
pub const __UINT16_FMTo__ = "ho";
pub const __UINT16_FMTu__ = "hu";
pub const __UINT16_FMTx__ = "hx";
pub const __UINT16_FMTX__ = "hX";
pub const __UINT16_C_SUFFIX__ = "";
pub const __UINT16_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65535, .decimal);
pub const __INT16_MAX__ = @as(c_int, 32767);
pub const __UINT32_TYPE__ = c_uint;
pub const __UINT32_FMTo__ = "o";
pub const __UINT32_FMTu__ = "u";
pub const __UINT32_FMTx__ = "x";
pub const __UINT32_FMTX__ = "X";
pub const __UINT32_C_SUFFIX__ = @compileError("unable to translate macro: undefined identifier `U`");
// (no file):224:9
pub const __UINT32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_uint, 4294967295, .decimal);
pub const __INT32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __UINT64_TYPE__ = c_ulong;
pub const __UINT64_FMTo__ = "lo";
pub const __UINT64_FMTu__ = "lu";
pub const __UINT64_FMTx__ = "lx";
pub const __UINT64_FMTX__ = "lX";
pub const __UINT64_C_SUFFIX__ = @compileError("unable to translate macro: undefined identifier `UL`");
// (no file):232:9
pub const __UINT64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __INT64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __INT_LEAST8_TYPE__ = i8;
pub const __INT_LEAST8_MAX__ = @as(c_int, 127);
pub const __INT_LEAST8_WIDTH__ = @as(c_int, 8);
pub const __INT_LEAST8_FMTd__ = "hhd";
pub const __INT_LEAST8_FMTi__ = "hhi";
pub const __UINT_LEAST8_TYPE__ = u8;
pub const __UINT_LEAST8_MAX__ = @as(c_int, 255);
pub const __UINT_LEAST8_FMTo__ = "hho";
pub const __UINT_LEAST8_FMTu__ = "hhu";
pub const __UINT_LEAST8_FMTx__ = "hhx";
pub const __UINT_LEAST8_FMTX__ = "hhX";
pub const __INT_LEAST16_TYPE__ = c_short;
pub const __INT_LEAST16_MAX__ = @as(c_int, 32767);
pub const __INT_LEAST16_WIDTH__ = @as(c_int, 16);
pub const __INT_LEAST16_FMTd__ = "hd";
pub const __INT_LEAST16_FMTi__ = "hi";
pub const __UINT_LEAST16_TYPE__ = c_ushort;
pub const __UINT_LEAST16_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65535, .decimal);
pub const __UINT_LEAST16_FMTo__ = "ho";
pub const __UINT_LEAST16_FMTu__ = "hu";
pub const __UINT_LEAST16_FMTx__ = "hx";
pub const __UINT_LEAST16_FMTX__ = "hX";
pub const __INT_LEAST32_TYPE__ = c_int;
pub const __INT_LEAST32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __INT_LEAST32_WIDTH__ = @as(c_int, 32);
pub const __INT_LEAST32_FMTd__ = "d";
pub const __INT_LEAST32_FMTi__ = "i";
pub const __UINT_LEAST32_TYPE__ = c_uint;
pub const __UINT_LEAST32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_uint, 4294967295, .decimal);
pub const __UINT_LEAST32_FMTo__ = "o";
pub const __UINT_LEAST32_FMTu__ = "u";
pub const __UINT_LEAST32_FMTx__ = "x";
pub const __UINT_LEAST32_FMTX__ = "X";
pub const __INT_LEAST64_TYPE__ = c_long;
pub const __INT_LEAST64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __INT_LEAST64_WIDTH__ = @as(c_int, 64);
pub const __INT_LEAST64_FMTd__ = "ld";
pub const __INT_LEAST64_FMTi__ = "li";
pub const __UINT_LEAST64_TYPE__ = c_ulong;
pub const __UINT_LEAST64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __UINT_LEAST64_FMTo__ = "lo";
pub const __UINT_LEAST64_FMTu__ = "lu";
pub const __UINT_LEAST64_FMTx__ = "lx";
pub const __UINT_LEAST64_FMTX__ = "lX";
pub const __INT_FAST8_TYPE__ = i8;
pub const __INT_FAST8_MAX__ = @as(c_int, 127);
pub const __INT_FAST8_WIDTH__ = @as(c_int, 8);
pub const __INT_FAST8_FMTd__ = "hhd";
pub const __INT_FAST8_FMTi__ = "hhi";
pub const __UINT_FAST8_TYPE__ = u8;
pub const __UINT_FAST8_MAX__ = @as(c_int, 255);
pub const __UINT_FAST8_FMTo__ = "hho";
pub const __UINT_FAST8_FMTu__ = "hhu";
pub const __UINT_FAST8_FMTx__ = "hhx";
pub const __UINT_FAST8_FMTX__ = "hhX";
pub const __INT_FAST16_TYPE__ = c_short;
pub const __INT_FAST16_MAX__ = @as(c_int, 32767);
pub const __INT_FAST16_WIDTH__ = @as(c_int, 16);
pub const __INT_FAST16_FMTd__ = "hd";
pub const __INT_FAST16_FMTi__ = "hi";
pub const __UINT_FAST16_TYPE__ = c_ushort;
pub const __UINT_FAST16_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65535, .decimal);
pub const __UINT_FAST16_FMTo__ = "ho";
pub const __UINT_FAST16_FMTu__ = "hu";
pub const __UINT_FAST16_FMTx__ = "hx";
pub const __UINT_FAST16_FMTX__ = "hX";
pub const __INT_FAST32_TYPE__ = c_int;
pub const __INT_FAST32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const __INT_FAST32_WIDTH__ = @as(c_int, 32);
pub const __INT_FAST32_FMTd__ = "d";
pub const __INT_FAST32_FMTi__ = "i";
pub const __UINT_FAST32_TYPE__ = c_uint;
pub const __UINT_FAST32_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_uint, 4294967295, .decimal);
pub const __UINT_FAST32_FMTo__ = "o";
pub const __UINT_FAST32_FMTu__ = "u";
pub const __UINT_FAST32_FMTx__ = "x";
pub const __UINT_FAST32_FMTX__ = "X";
pub const __INT_FAST64_TYPE__ = c_long;
pub const __INT_FAST64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_long, 9223372036854775807, .decimal);
pub const __INT_FAST64_WIDTH__ = @as(c_int, 64);
pub const __INT_FAST64_FMTd__ = "ld";
pub const __INT_FAST64_FMTi__ = "li";
pub const __UINT_FAST64_TYPE__ = c_ulong;
pub const __UINT_FAST64_MAX__ = @import("std").zig.c_translation.promoteIntLiteral(c_ulong, 18446744073709551615, .decimal);
pub const __UINT_FAST64_FMTo__ = "lo";
pub const __UINT_FAST64_FMTu__ = "lu";
pub const __UINT_FAST64_FMTx__ = "lx";
pub const __UINT_FAST64_FMTX__ = "lX";
pub const __USER_LABEL_PREFIX__ = "";
pub const __FINITE_MATH_ONLY__ = @as(c_int, 0);
pub const __GNUC_STDC_INLINE__ = @as(c_int, 1);
pub const __GCC_ATOMIC_TEST_AND_SET_TRUEVAL = @as(c_int, 1);
pub const __GCC_DESTRUCTIVE_SIZE = @as(c_int, 64);
pub const __GCC_CONSTRUCTIVE_SIZE = @as(c_int, 64);
pub const __CLANG_ATOMIC_BOOL_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_CHAR_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_CHAR16_T_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_CHAR32_T_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_WCHAR_T_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_SHORT_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_INT_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_LONG_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_LLONG_LOCK_FREE = @as(c_int, 2);
pub const __CLANG_ATOMIC_POINTER_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_BOOL_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_CHAR_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_CHAR16_T_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_CHAR32_T_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_WCHAR_T_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_SHORT_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_INT_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_LONG_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_LLONG_LOCK_FREE = @as(c_int, 2);
pub const __GCC_ATOMIC_POINTER_LOCK_FREE = @as(c_int, 2);
pub const __NO_INLINE__ = @as(c_int, 1);
pub const __PIC__ = @as(c_int, 2);
pub const __pic__ = @as(c_int, 2);
pub const __FLT_RADIX__ = @as(c_int, 2);
pub const __DECIMAL_DIG__ = __LDBL_DECIMAL_DIG__;
pub const __SSP_STRONG__ = @as(c_int, 2);
pub const __ELF__ = @as(c_int, 1);
pub const __GCC_ASM_FLAG_OUTPUTS__ = @as(c_int, 1);
pub const __code_model_small__ = @as(c_int, 1);
pub const __amd64__ = @as(c_int, 1);
pub const __amd64 = @as(c_int, 1);
pub const __x86_64 = @as(c_int, 1);
pub const __x86_64__ = @as(c_int, 1);
pub const __SEG_GS = @as(c_int, 1);
pub const __SEG_FS = @as(c_int, 1);
pub const __seg_gs = @compileError("unable to translate macro: undefined identifier `address_space`");
// (no file):366:9
pub const __seg_fs = @compileError("unable to translate macro: undefined identifier `address_space`");
// (no file):367:9
pub const __znver4 = @as(c_int, 1);
pub const __znver4__ = @as(c_int, 1);
pub const __tune_znver4__ = @as(c_int, 1);
pub const __REGISTER_PREFIX__ = "";
pub const __NO_MATH_INLINES = @as(c_int, 1);
pub const __AES__ = @as(c_int, 1);
pub const __VAES__ = @as(c_int, 1);
pub const __PCLMUL__ = @as(c_int, 1);
pub const __VPCLMULQDQ__ = @as(c_int, 1);
pub const __LAHF_SAHF__ = @as(c_int, 1);
pub const __LZCNT__ = @as(c_int, 1);
pub const __RDRND__ = @as(c_int, 1);
pub const __FSGSBASE__ = @as(c_int, 1);
pub const __BMI__ = @as(c_int, 1);
pub const __BMI2__ = @as(c_int, 1);
pub const __POPCNT__ = @as(c_int, 1);
pub const __PRFCHW__ = @as(c_int, 1);
pub const __RDSEED__ = @as(c_int, 1);
pub const __ADX__ = @as(c_int, 1);
pub const __MWAITX__ = @as(c_int, 1);
pub const __MOVBE__ = @as(c_int, 1);
pub const __SSE4A__ = @as(c_int, 1);
pub const __FMA__ = @as(c_int, 1);
pub const __F16C__ = @as(c_int, 1);
pub const __GFNI__ = @as(c_int, 1);
pub const __EVEX512__ = @as(c_int, 1);
pub const __AVX512CD__ = @as(c_int, 1);
pub const __AVX512VPOPCNTDQ__ = @as(c_int, 1);
pub const __AVX512VNNI__ = @as(c_int, 1);
pub const __AVX512BF16__ = @as(c_int, 1);
pub const __AVX512DQ__ = @as(c_int, 1);
pub const __AVX512BITALG__ = @as(c_int, 1);
pub const __AVX512BW__ = @as(c_int, 1);
pub const __AVX512VL__ = @as(c_int, 1);
pub const __EVEX256__ = @as(c_int, 1);
pub const __AVX512VBMI__ = @as(c_int, 1);
pub const __AVX512VBMI2__ = @as(c_int, 1);
pub const __AVX512IFMA__ = @as(c_int, 1);
pub const __SHA__ = @as(c_int, 1);
pub const __FXSR__ = @as(c_int, 1);
pub const __XSAVE__ = @as(c_int, 1);
pub const __XSAVEOPT__ = @as(c_int, 1);
pub const __XSAVEC__ = @as(c_int, 1);
pub const __XSAVES__ = @as(c_int, 1);
pub const __PKU__ = @as(c_int, 1);
pub const __CLFLUSHOPT__ = @as(c_int, 1);
pub const __CLWB__ = @as(c_int, 1);
pub const __WBNOINVD__ = @as(c_int, 1);
pub const __SHSTK__ = @as(c_int, 1);
pub const __CLZERO__ = @as(c_int, 1);
pub const __RDPID__ = @as(c_int, 1);
pub const __RDPRU__ = @as(c_int, 1);
pub const __INVPCID__ = @as(c_int, 1);
pub const __CRC32__ = @as(c_int, 1);
pub const __AVX512F__ = @as(c_int, 1);
pub const __AVX2__ = @as(c_int, 1);
pub const __AVX__ = @as(c_int, 1);
pub const __SSE4_2__ = @as(c_int, 1);
pub const __SSE4_1__ = @as(c_int, 1);
pub const __SSSE3__ = @as(c_int, 1);
pub const __SSE3__ = @as(c_int, 1);
pub const __SSE2__ = @as(c_int, 1);
pub const __SSE2_MATH__ = @as(c_int, 1);
pub const __SSE__ = @as(c_int, 1);
pub const __SSE_MATH__ = @as(c_int, 1);
pub const __MMX__ = @as(c_int, 1);
pub const __GCC_HAVE_SYNC_COMPARE_AND_SWAP_1 = @as(c_int, 1);
pub const __GCC_HAVE_SYNC_COMPARE_AND_SWAP_2 = @as(c_int, 1);
pub const __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4 = @as(c_int, 1);
pub const __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8 = @as(c_int, 1);
pub const __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16 = @as(c_int, 1);
pub const __SIZEOF_FLOAT128__ = @as(c_int, 16);
pub const unix = @as(c_int, 1);
pub const __unix = @as(c_int, 1);
pub const __unix__ = @as(c_int, 1);
pub const linux = @as(c_int, 1);
pub const __linux = @as(c_int, 1);
pub const __linux__ = @as(c_int, 1);
pub const __gnu_linux__ = @as(c_int, 1);
pub const __FLOAT128__ = @as(c_int, 1);
pub const __STDC__ = @as(c_int, 1);
pub const __STDC_HOSTED__ = @as(c_int, 1);
pub const __STDC_VERSION__ = @as(c_long, 201710);
pub const __STDC_UTF_16__ = @as(c_int, 1);
pub const __STDC_UTF_32__ = @as(c_int, 1);
pub const __STDC_EMBED_NOT_FOUND__ = @as(c_int, 0);
pub const __STDC_EMBED_FOUND__ = @as(c_int, 1);
pub const __STDC_EMBED_EMPTY__ = @as(c_int, 2);
pub const __GLIBC_MINOR__ = @as(c_int, 41);
pub const _DEBUG = @as(c_int, 1);
pub const __GCC_HAVE_DWARF2_CFI_ASM = @as(c_int, 1);
pub const __GSL_LINALG_H__ = "";
pub const __GLIBC_INTERNAL_STARTING_HEADER_IMPLEMENTATION = "";
pub const _FEATURES_H = @as(c_int, 1);
pub const __KERNEL_STRICT_NAMES = "";
pub inline fn __GNUC_PREREQ(maj: anytype, min: anytype) @TypeOf(((__GNUC__ << @as(c_int, 16)) + __GNUC_MINOR__) >= ((maj << @as(c_int, 16)) + min)) {
    _ = &maj;
    _ = &min;
    return ((__GNUC__ << @as(c_int, 16)) + __GNUC_MINOR__) >= ((maj << @as(c_int, 16)) + min);
}
pub inline fn __glibc_clang_prereq(maj: anytype, min: anytype) @TypeOf(((__clang_major__ << @as(c_int, 16)) + __clang_minor__) >= ((maj << @as(c_int, 16)) + min)) {
    _ = &maj;
    _ = &min;
    return ((__clang_major__ << @as(c_int, 16)) + __clang_minor__) >= ((maj << @as(c_int, 16)) + min);
}
pub const __GLIBC_USE = @compileError("unable to translate macro: undefined identifier `__GLIBC_USE_`");
// /usr/include/features.h:191:9
pub const _DEFAULT_SOURCE = @as(c_int, 1);
pub const __GLIBC_USE_ISOC2Y = @as(c_int, 0);
pub const __GLIBC_USE_ISOC23 = @as(c_int, 0);
pub const __USE_ISOC11 = @as(c_int, 1);
pub const __USE_ISOC99 = @as(c_int, 1);
pub const __USE_ISOC95 = @as(c_int, 1);
pub const __USE_POSIX_IMPLICITLY = @as(c_int, 1);
pub const _POSIX_SOURCE = @as(c_int, 1);
pub const _POSIX_C_SOURCE = @as(c_long, 200809);
pub const __USE_POSIX = @as(c_int, 1);
pub const __USE_POSIX2 = @as(c_int, 1);
pub const __USE_POSIX199309 = @as(c_int, 1);
pub const __USE_POSIX199506 = @as(c_int, 1);
pub const __USE_XOPEN2K = @as(c_int, 1);
pub const __USE_XOPEN2K8 = @as(c_int, 1);
pub const _ATFILE_SOURCE = @as(c_int, 1);
pub const __WORDSIZE = @as(c_int, 64);
pub const __WORDSIZE_TIME64_COMPAT32 = @as(c_int, 1);
pub const __SYSCALL_WORDSIZE = @as(c_int, 64);
pub const __TIMESIZE = __WORDSIZE;
pub const __USE_TIME_BITS64 = @as(c_int, 1);
pub const __USE_MISC = @as(c_int, 1);
pub const __USE_ATFILE = @as(c_int, 1);
pub const __USE_FORTIFY_LEVEL = @as(c_int, 0);
pub const __GLIBC_USE_DEPRECATED_GETS = @as(c_int, 0);
pub const __GLIBC_USE_DEPRECATED_SCANF = @as(c_int, 0);
pub const __GLIBC_USE_C23_STRTOL = @as(c_int, 0);
pub const _STDC_PREDEF_H = @as(c_int, 1);
pub const __STDC_IEC_559__ = @as(c_int, 1);
pub const __STDC_IEC_60559_BFP__ = @as(c_long, 201404);
pub const __STDC_IEC_559_COMPLEX__ = @as(c_int, 1);
pub const __STDC_IEC_60559_COMPLEX__ = @as(c_long, 201404);
pub const __STDC_ISO_10646__ = @as(c_long, 201706);
pub const __GNU_LIBRARY__ = @as(c_int, 6);
pub const __GLIBC__ = @as(c_int, 2);
pub inline fn __GLIBC_PREREQ(maj: anytype, min: anytype) @TypeOf(((__GLIBC__ << @as(c_int, 16)) + __GLIBC_MINOR__) >= ((maj << @as(c_int, 16)) + min)) {
    _ = &maj;
    _ = &min;
    return ((__GLIBC__ << @as(c_int, 16)) + __GLIBC_MINOR__) >= ((maj << @as(c_int, 16)) + min);
}
pub const _SYS_CDEFS_H = @as(c_int, 1);
pub const __glibc_has_attribute = @compileError("unable to translate macro: undefined identifier `__has_attribute`");
// /usr/include/sys/cdefs.h:45:10
pub inline fn __glibc_has_builtin(name: anytype) @TypeOf(__has_builtin(name)) {
    _ = &name;
    return __has_builtin(name);
}
pub const __glibc_has_extension = @compileError("unable to translate macro: undefined identifier `__has_extension`");
// /usr/include/sys/cdefs.h:55:10
pub const __LEAF = "";
pub const __LEAF_ATTR = "";
pub const __THROW = @compileError("unable to translate macro: undefined identifier `__nothrow__`");
// /usr/include/sys/cdefs.h:79:11
pub const __THROWNL = @compileError("unable to translate macro: undefined identifier `__nothrow__`");
// /usr/include/sys/cdefs.h:80:11
pub const __NTH = @compileError("unable to translate macro: undefined identifier `__nothrow__`");
// /usr/include/sys/cdefs.h:81:11
pub const __NTHNL = @compileError("unable to translate macro: undefined identifier `__nothrow__`");
// /usr/include/sys/cdefs.h:82:11
pub const __COLD = @compileError("unable to translate macro: undefined identifier `__cold__`");
// /usr/include/sys/cdefs.h:102:11
pub inline fn __P(args: anytype) @TypeOf(args) {
    _ = &args;
    return args;
}
pub inline fn __PMT(args: anytype) @TypeOf(args) {
    _ = &args;
    return args;
}
pub const __CONCAT = @compileError("unable to translate C expr: unexpected token '##'");
// /usr/include/sys/cdefs.h:131:9
pub const __STRING = @compileError("unable to translate C expr: unexpected token '#'");
// /usr/include/sys/cdefs.h:132:9
pub const __ptr_t = ?*anyopaque;
pub const __BEGIN_DECLS = "";
pub const __END_DECLS = "";
pub const __attribute_overloadable__ = @compileError("unable to translate macro: undefined identifier `__overloadable__`");
// /usr/include/sys/cdefs.h:151:10
pub inline fn __bos(ptr: anytype) @TypeOf(__builtin_object_size(ptr, __USE_FORTIFY_LEVEL > @as(c_int, 1))) {
    _ = &ptr;
    return __builtin_object_size(ptr, __USE_FORTIFY_LEVEL > @as(c_int, 1));
}
pub inline fn __bos0(ptr: anytype) @TypeOf(__builtin_object_size(ptr, @as(c_int, 0))) {
    _ = &ptr;
    return __builtin_object_size(ptr, @as(c_int, 0));
}
pub inline fn __glibc_objsize0(__o: anytype) @TypeOf(__bos0(__o)) {
    _ = &__o;
    return __bos0(__o);
}
pub inline fn __glibc_objsize(__o: anytype) @TypeOf(__bos(__o)) {
    _ = &__o;
    return __bos(__o);
}
pub const __warnattr = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:370:10
pub const __errordecl = @compileError("unable to translate C expr: unexpected token 'extern'");
// /usr/include/sys/cdefs.h:371:10
pub const __flexarr = @compileError("unable to translate C expr: unexpected token '['");
// /usr/include/sys/cdefs.h:379:10
pub const __glibc_c99_flexarr_available = @as(c_int, 1);
pub const __REDIRECT = @compileError("unable to translate C expr: unexpected token '__asm__'");
// /usr/include/sys/cdefs.h:410:10
pub const __REDIRECT_NTH = @compileError("unable to translate C expr: unexpected token '__asm__'");
// /usr/include/sys/cdefs.h:417:11
pub const __REDIRECT_NTHNL = @compileError("unable to translate C expr: unexpected token '__asm__'");
// /usr/include/sys/cdefs.h:419:11
pub const __ASMNAME = @compileError("unable to translate C expr: unexpected token ','");
// /usr/include/sys/cdefs.h:422:10
pub inline fn __ASMNAME2(prefix: anytype, cname: anytype) @TypeOf(__STRING(prefix) ++ cname) {
    _ = &prefix;
    _ = &cname;
    return __STRING(prefix) ++ cname;
}
pub const __REDIRECT_FORTIFY = __REDIRECT;
pub const __REDIRECT_FORTIFY_NTH = __REDIRECT_NTH;
pub const __attribute_malloc__ = @compileError("unable to translate macro: undefined identifier `__malloc__`");
// /usr/include/sys/cdefs.h:452:10
pub const __attribute_alloc_size__ = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:463:10
pub const __attribute_alloc_align__ = @compileError("unable to translate macro: undefined identifier `__alloc_align__`");
// /usr/include/sys/cdefs.h:469:10
pub const __attribute_pure__ = @compileError("unable to translate macro: undefined identifier `__pure__`");
// /usr/include/sys/cdefs.h:479:10
pub const __attribute_const__ = @compileError("unable to translate C expr: unexpected token '__attribute__'");
// /usr/include/sys/cdefs.h:486:10
pub const __attribute_maybe_unused__ = @compileError("unable to translate macro: undefined identifier `__unused__`");
// /usr/include/sys/cdefs.h:492:10
pub const __attribute_used__ = @compileError("unable to translate macro: undefined identifier `__used__`");
// /usr/include/sys/cdefs.h:501:10
pub const __attribute_noinline__ = @compileError("unable to translate macro: undefined identifier `__noinline__`");
// /usr/include/sys/cdefs.h:502:10
pub const __attribute_deprecated__ = @compileError("unable to translate macro: undefined identifier `__deprecated__`");
// /usr/include/sys/cdefs.h:510:10
pub const __attribute_deprecated_msg__ = @compileError("unable to translate macro: undefined identifier `__deprecated__`");
// /usr/include/sys/cdefs.h:520:10
pub const __attribute_format_arg__ = @compileError("unable to translate macro: undefined identifier `__format_arg__`");
// /usr/include/sys/cdefs.h:533:10
pub const __attribute_format_strfmon__ = @compileError("unable to translate macro: undefined identifier `__format__`");
// /usr/include/sys/cdefs.h:543:10
pub const __attribute_nonnull__ = @compileError("unable to translate macro: undefined identifier `__nonnull__`");
// /usr/include/sys/cdefs.h:555:11
pub inline fn __nonnull(params: anytype) @TypeOf(__attribute_nonnull__(params)) {
    _ = &params;
    return __attribute_nonnull__(params);
}
pub const __returns_nonnull = @compileError("unable to translate macro: undefined identifier `__returns_nonnull__`");
// /usr/include/sys/cdefs.h:568:10
pub const __attribute_warn_unused_result__ = @compileError("unable to translate macro: undefined identifier `__warn_unused_result__`");
// /usr/include/sys/cdefs.h:577:10
pub const __wur = "";
pub const __always_inline = @compileError("unable to translate macro: undefined identifier `__always_inline__`");
// /usr/include/sys/cdefs.h:595:10
pub const __attribute_artificial__ = @compileError("unable to translate macro: undefined identifier `__artificial__`");
// /usr/include/sys/cdefs.h:604:10
pub const __extern_inline = @compileError("unable to translate macro: undefined identifier `__gnu_inline__`");
// /usr/include/sys/cdefs.h:622:11
pub const __extern_always_inline = @compileError("unable to translate macro: undefined identifier `__gnu_inline__`");
// /usr/include/sys/cdefs.h:623:11
pub const __fortify_function = __extern_always_inline ++ __attribute_artificial__;
pub const __restrict_arr = @compileError("unable to translate C expr: unexpected token '__restrict'");
// /usr/include/sys/cdefs.h:666:10
pub inline fn __glibc_unlikely(cond: anytype) @TypeOf(__builtin_expect(cond, @as(c_int, 0))) {
    _ = &cond;
    return __builtin_expect(cond, @as(c_int, 0));
}
pub inline fn __glibc_likely(cond: anytype) @TypeOf(__builtin_expect(cond, @as(c_int, 1))) {
    _ = &cond;
    return __builtin_expect(cond, @as(c_int, 1));
}
pub const __attribute_nonstring__ = "";
pub const __attribute_copy__ = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:715:10
pub const __LDOUBLE_REDIRECTS_TO_FLOAT128_ABI = @as(c_int, 0);
pub inline fn __LDBL_REDIR1(name: anytype, proto: anytype, alias: anytype) @TypeOf(name ++ proto) {
    _ = &name;
    _ = &proto;
    _ = &alias;
    return name ++ proto;
}
pub inline fn __LDBL_REDIR(name: anytype, proto: anytype) @TypeOf(name ++ proto) {
    _ = &name;
    _ = &proto;
    return name ++ proto;
}
pub inline fn __LDBL_REDIR1_NTH(name: anytype, proto: anytype, alias: anytype) @TypeOf(name ++ proto ++ __THROW) {
    _ = &name;
    _ = &proto;
    _ = &alias;
    return name ++ proto ++ __THROW;
}
pub inline fn __LDBL_REDIR_NTH(name: anytype, proto: anytype) @TypeOf(name ++ proto ++ __THROW) {
    _ = &name;
    _ = &proto;
    return name ++ proto ++ __THROW;
}
pub const __LDBL_REDIR2_DECL = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:792:10
pub const __LDBL_REDIR_DECL = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:793:10
pub inline fn __REDIRECT_LDBL(name: anytype, proto: anytype, alias: anytype) @TypeOf(__REDIRECT(name, proto, alias)) {
    _ = &name;
    _ = &proto;
    _ = &alias;
    return __REDIRECT(name, proto, alias);
}
pub inline fn __REDIRECT_NTH_LDBL(name: anytype, proto: anytype, alias: anytype) @TypeOf(__REDIRECT_NTH(name, proto, alias)) {
    _ = &name;
    _ = &proto;
    _ = &alias;
    return __REDIRECT_NTH(name, proto, alias);
}
pub const __glibc_macro_warning1 = @compileError("unable to translate macro: undefined identifier `_Pragma`");
// /usr/include/sys/cdefs.h:807:10
pub const __glibc_macro_warning = @compileError("unable to translate macro: undefined identifier `GCC`");
// /usr/include/sys/cdefs.h:808:10
pub const __HAVE_GENERIC_SELECTION = @as(c_int, 1);
pub const __fortified_attr_access = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:853:11
pub const __attr_access = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:854:11
pub const __attr_access_none = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:855:11
pub const __attr_dealloc = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/sys/cdefs.h:865:10
pub const __attr_dealloc_free = "";
pub const __attribute_returns_twice__ = @compileError("unable to translate macro: undefined identifier `__returns_twice__`");
// /usr/include/sys/cdefs.h:872:10
pub const __attribute_struct_may_alias__ = @compileError("unable to translate macro: undefined identifier `__may_alias__`");
// /usr/include/sys/cdefs.h:881:10
pub const __stub___compat_bdflush = "";
pub const __stub_chflags = "";
pub const __stub_fchflags = "";
pub const __stub_gtty = "";
pub const __stub_revoke = "";
pub const __stub_setlogin = "";
pub const __stub_sigreturn = "";
pub const __stub_stty = "";
pub const __GLIBC_USE_LIB_EXT2 = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_BFP_EXT = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_BFP_EXT_C23 = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_EXT = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_FUNCS_EXT = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_FUNCS_EXT_C23 = @as(c_int, 0);
pub const __GLIBC_USE_IEC_60559_TYPES_EXT = @as(c_int, 0);
pub const __need_size_t = "";
pub const __need_wchar_t = "";
pub const __need_NULL = "";
pub const _SIZE_T = "";
pub const _WCHAR_T = "";
pub const NULL = @import("std").zig.c_translation.cast(?*anyopaque, @as(c_int, 0));
pub const _STDLIB_H = @as(c_int, 1);
pub const WNOHANG = @as(c_int, 1);
pub const WUNTRACED = @as(c_int, 2);
pub const WSTOPPED = @as(c_int, 2);
pub const WEXITED = @as(c_int, 4);
pub const WCONTINUED = @as(c_int, 8);
pub const WNOWAIT = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0x01000000, .hex);
pub const __WNOTHREAD = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0x20000000, .hex);
pub const __WALL = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0x40000000, .hex);
pub const __WCLONE = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0x80000000, .hex);
pub inline fn __WEXITSTATUS(status: anytype) @TypeOf((status & @import("std").zig.c_translation.promoteIntLiteral(c_int, 0xff00, .hex)) >> @as(c_int, 8)) {
    _ = &status;
    return (status & @import("std").zig.c_translation.promoteIntLiteral(c_int, 0xff00, .hex)) >> @as(c_int, 8);
}
pub inline fn __WTERMSIG(status: anytype) @TypeOf(status & @as(c_int, 0x7f)) {
    _ = &status;
    return status & @as(c_int, 0x7f);
}
pub inline fn __WSTOPSIG(status: anytype) @TypeOf(__WEXITSTATUS(status)) {
    _ = &status;
    return __WEXITSTATUS(status);
}
pub inline fn __WIFEXITED(status: anytype) @TypeOf(__WTERMSIG(status) == @as(c_int, 0)) {
    _ = &status;
    return __WTERMSIG(status) == @as(c_int, 0);
}
pub inline fn __WIFSIGNALED(status: anytype) @TypeOf((@import("std").zig.c_translation.cast(i8, (status & @as(c_int, 0x7f)) + @as(c_int, 1)) >> @as(c_int, 1)) > @as(c_int, 0)) {
    _ = &status;
    return (@import("std").zig.c_translation.cast(i8, (status & @as(c_int, 0x7f)) + @as(c_int, 1)) >> @as(c_int, 1)) > @as(c_int, 0);
}
pub inline fn __WIFSTOPPED(status: anytype) @TypeOf((status & @as(c_int, 0xff)) == @as(c_int, 0x7f)) {
    _ = &status;
    return (status & @as(c_int, 0xff)) == @as(c_int, 0x7f);
}
pub inline fn __WIFCONTINUED(status: anytype) @TypeOf(status == __W_CONTINUED) {
    _ = &status;
    return status == __W_CONTINUED;
}
pub inline fn __WCOREDUMP(status: anytype) @TypeOf(status & __WCOREFLAG) {
    _ = &status;
    return status & __WCOREFLAG;
}
pub inline fn __W_EXITCODE(ret: anytype, sig: anytype) @TypeOf((ret << @as(c_int, 8)) | sig) {
    _ = &ret;
    _ = &sig;
    return (ret << @as(c_int, 8)) | sig;
}
pub inline fn __W_STOPCODE(sig: anytype) @TypeOf((sig << @as(c_int, 8)) | @as(c_int, 0x7f)) {
    _ = &sig;
    return (sig << @as(c_int, 8)) | @as(c_int, 0x7f);
}
pub const __W_CONTINUED = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0xffff, .hex);
pub const __WCOREFLAG = @as(c_int, 0x80);
pub inline fn WEXITSTATUS(status: anytype) @TypeOf(__WEXITSTATUS(status)) {
    _ = &status;
    return __WEXITSTATUS(status);
}
pub inline fn WTERMSIG(status: anytype) @TypeOf(__WTERMSIG(status)) {
    _ = &status;
    return __WTERMSIG(status);
}
pub inline fn WSTOPSIG(status: anytype) @TypeOf(__WSTOPSIG(status)) {
    _ = &status;
    return __WSTOPSIG(status);
}
pub inline fn WIFEXITED(status: anytype) @TypeOf(__WIFEXITED(status)) {
    _ = &status;
    return __WIFEXITED(status);
}
pub inline fn WIFSIGNALED(status: anytype) @TypeOf(__WIFSIGNALED(status)) {
    _ = &status;
    return __WIFSIGNALED(status);
}
pub inline fn WIFSTOPPED(status: anytype) @TypeOf(__WIFSTOPPED(status)) {
    _ = &status;
    return __WIFSTOPPED(status);
}
pub inline fn WIFCONTINUED(status: anytype) @TypeOf(__WIFCONTINUED(status)) {
    _ = &status;
    return __WIFCONTINUED(status);
}
pub const _BITS_FLOATN_H = "";
pub const __HAVE_FLOAT128 = @as(c_int, 1);
pub const __HAVE_DISTINCT_FLOAT128 = @as(c_int, 1);
pub const __HAVE_FLOAT64X = @as(c_int, 1);
pub const __HAVE_FLOAT64X_LONG_DOUBLE = @as(c_int, 1);
pub const __f128 = @compileError("unable to translate macro: undefined identifier `q`");
// /usr/include/bits/floatn.h:70:12
pub const __CFLOAT128 = __cfloat128;
pub const __builtin_signbitf128 = __signbitf128;
pub const _BITS_FLOATN_COMMON_H = "";
pub const __HAVE_FLOAT16 = @as(c_int, 0);
pub const __HAVE_FLOAT32 = @as(c_int, 1);
pub const __HAVE_FLOAT64 = @as(c_int, 1);
pub const __HAVE_FLOAT32X = @as(c_int, 1);
pub const __HAVE_FLOAT128X = @as(c_int, 0);
pub const __HAVE_DISTINCT_FLOAT16 = __HAVE_FLOAT16;
pub const __HAVE_DISTINCT_FLOAT32 = @as(c_int, 0);
pub const __HAVE_DISTINCT_FLOAT64 = @as(c_int, 0);
pub const __HAVE_DISTINCT_FLOAT32X = @as(c_int, 0);
pub const __HAVE_DISTINCT_FLOAT64X = @as(c_int, 0);
pub const __HAVE_DISTINCT_FLOAT128X = __HAVE_FLOAT128X;
pub const __HAVE_FLOAT128_UNLIKE_LDBL = (__HAVE_DISTINCT_FLOAT128 != 0) and (__LDBL_MANT_DIG__ != @as(c_int, 113));
pub const __HAVE_FLOATN_NOT_TYPEDEF = @as(c_int, 0);
pub const __f32 = @import("std").zig.c_translation.Macros.F_SUFFIX;
pub inline fn __f64(x: anytype) @TypeOf(x) {
    _ = &x;
    return x;
}
pub inline fn __f32x(x: anytype) @TypeOf(x) {
    _ = &x;
    return x;
}
pub const __f64x = @import("std").zig.c_translation.Macros.L_SUFFIX;
pub const __CFLOAT32 = @compileError("unable to translate: TODO _Complex");
// /usr/include/bits/floatn-common.h:149:12
pub const __CFLOAT64 = @compileError("unable to translate: TODO _Complex");
// /usr/include/bits/floatn-common.h:160:13
pub const __CFLOAT32X = @compileError("unable to translate: TODO _Complex");
// /usr/include/bits/floatn-common.h:169:12
pub const __CFLOAT64X = @compileError("unable to translate: TODO _Complex");
// /usr/include/bits/floatn-common.h:178:13
pub inline fn __builtin_huge_valf32() @TypeOf(__builtin_huge_valf()) {
    return __builtin_huge_valf();
}
pub inline fn __builtin_inff32() @TypeOf(__builtin_inff()) {
    return __builtin_inff();
}
pub inline fn __builtin_nanf32(x: anytype) @TypeOf(__builtin_nanf(x)) {
    _ = &x;
    return __builtin_nanf(x);
}
pub const __builtin_nansf32 = @compileError("unable to translate macro: undefined identifier `__builtin_nansf`");
// /usr/include/bits/floatn-common.h:221:12
pub const __builtin_huge_valf64 = @compileError("unable to translate macro: undefined identifier `__builtin_huge_val`");
// /usr/include/bits/floatn-common.h:255:13
pub const __builtin_inff64 = @compileError("unable to translate macro: undefined identifier `__builtin_inf`");
// /usr/include/bits/floatn-common.h:256:13
pub const __builtin_nanf64 = @compileError("unable to translate macro: undefined identifier `__builtin_nan`");
// /usr/include/bits/floatn-common.h:257:13
pub const __builtin_nansf64 = @compileError("unable to translate macro: undefined identifier `__builtin_nans`");
// /usr/include/bits/floatn-common.h:258:13
pub const __builtin_huge_valf32x = @compileError("unable to translate macro: undefined identifier `__builtin_huge_val`");
// /usr/include/bits/floatn-common.h:272:12
pub const __builtin_inff32x = @compileError("unable to translate macro: undefined identifier `__builtin_inf`");
// /usr/include/bits/floatn-common.h:273:12
pub const __builtin_nanf32x = @compileError("unable to translate macro: undefined identifier `__builtin_nan`");
// /usr/include/bits/floatn-common.h:274:12
pub const __builtin_nansf32x = @compileError("unable to translate macro: undefined identifier `__builtin_nans`");
// /usr/include/bits/floatn-common.h:275:12
pub const __builtin_huge_valf64x = @compileError("unable to translate macro: undefined identifier `__builtin_huge_vall`");
// /usr/include/bits/floatn-common.h:289:13
pub const __builtin_inff64x = @compileError("unable to translate macro: undefined identifier `__builtin_infl`");
// /usr/include/bits/floatn-common.h:290:13
pub const __builtin_nanf64x = @compileError("unable to translate macro: undefined identifier `__builtin_nanl`");
// /usr/include/bits/floatn-common.h:291:13
pub const __builtin_nansf64x = @compileError("unable to translate macro: undefined identifier `__builtin_nansl`");
// /usr/include/bits/floatn-common.h:292:13
pub const __ldiv_t_defined = @as(c_int, 1);
pub const __lldiv_t_defined = @as(c_int, 1);
pub const RAND_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const EXIT_FAILURE = @as(c_int, 1);
pub const EXIT_SUCCESS = @as(c_int, 0);
pub const MB_CUR_MAX = __ctype_get_mb_cur_max();
pub const _SYS_TYPES_H = @as(c_int, 1);
pub const _BITS_TYPES_H = @as(c_int, 1);
pub const __S16_TYPE = c_short;
pub const __U16_TYPE = c_ushort;
pub const __S32_TYPE = c_int;
pub const __U32_TYPE = c_uint;
pub const __SLONGWORD_TYPE = c_long;
pub const __ULONGWORD_TYPE = c_ulong;
pub const __SQUAD_TYPE = c_long;
pub const __UQUAD_TYPE = c_ulong;
pub const __SWORD_TYPE = c_long;
pub const __UWORD_TYPE = c_ulong;
pub const __SLONG32_TYPE = c_int;
pub const __ULONG32_TYPE = c_uint;
pub const __S64_TYPE = c_long;
pub const __U64_TYPE = c_ulong;
pub const __STD_TYPE = @compileError("unable to translate C expr: unexpected token 'typedef'");
// /usr/include/bits/types.h:137:10
pub const _BITS_TYPESIZES_H = @as(c_int, 1);
pub const __SYSCALL_SLONG_TYPE = __SLONGWORD_TYPE;
pub const __SYSCALL_ULONG_TYPE = __ULONGWORD_TYPE;
pub const __DEV_T_TYPE = __UQUAD_TYPE;
pub const __UID_T_TYPE = __U32_TYPE;
pub const __GID_T_TYPE = __U32_TYPE;
pub const __INO_T_TYPE = __SYSCALL_ULONG_TYPE;
pub const __INO64_T_TYPE = __UQUAD_TYPE;
pub const __MODE_T_TYPE = __U32_TYPE;
pub const __NLINK_T_TYPE = __SYSCALL_ULONG_TYPE;
pub const __FSWORD_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __OFF_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __OFF64_T_TYPE = __SQUAD_TYPE;
pub const __PID_T_TYPE = __S32_TYPE;
pub const __RLIM_T_TYPE = __SYSCALL_ULONG_TYPE;
pub const __RLIM64_T_TYPE = __UQUAD_TYPE;
pub const __BLKCNT_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __BLKCNT64_T_TYPE = __SQUAD_TYPE;
pub const __FSBLKCNT_T_TYPE = __SYSCALL_ULONG_TYPE;
pub const __FSBLKCNT64_T_TYPE = __UQUAD_TYPE;
pub const __FSFILCNT_T_TYPE = __SYSCALL_ULONG_TYPE;
pub const __FSFILCNT64_T_TYPE = __UQUAD_TYPE;
pub const __ID_T_TYPE = __U32_TYPE;
pub const __CLOCK_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __TIME_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __USECONDS_T_TYPE = __U32_TYPE;
pub const __SUSECONDS_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __SUSECONDS64_T_TYPE = __SQUAD_TYPE;
pub const __DADDR_T_TYPE = __S32_TYPE;
pub const __KEY_T_TYPE = __S32_TYPE;
pub const __CLOCKID_T_TYPE = __S32_TYPE;
pub const __TIMER_T_TYPE = ?*anyopaque;
pub const __BLKSIZE_T_TYPE = __SYSCALL_SLONG_TYPE;
pub const __FSID_T_TYPE = @compileError("unable to translate macro: undefined identifier `__val`");
// /usr/include/bits/typesizes.h:73:9
pub const __SSIZE_T_TYPE = __SWORD_TYPE;
pub const __CPU_MASK_TYPE = __SYSCALL_ULONG_TYPE;
pub const __OFF_T_MATCHES_OFF64_T = @as(c_int, 1);
pub const __INO_T_MATCHES_INO64_T = @as(c_int, 1);
pub const __RLIM_T_MATCHES_RLIM64_T = @as(c_int, 1);
pub const __STATFS_MATCHES_STATFS64 = @as(c_int, 1);
pub const __KERNEL_OLD_TIMEVAL_MATCHES_TIMEVAL64 = @as(c_int, 1);
pub const __FD_SETSIZE = @as(c_int, 1024);
pub const _BITS_TIME64_H = @as(c_int, 1);
pub const __TIME64_T_TYPE = __TIME_T_TYPE;
pub const __u_char_defined = "";
pub const __ino_t_defined = "";
pub const __dev_t_defined = "";
pub const __gid_t_defined = "";
pub const __mode_t_defined = "";
pub const __nlink_t_defined = "";
pub const __uid_t_defined = "";
pub const __off_t_defined = "";
pub const __pid_t_defined = "";
pub const __id_t_defined = "";
pub const __ssize_t_defined = "";
pub const __daddr_t_defined = "";
pub const __key_t_defined = "";
pub const __clock_t_defined = @as(c_int, 1);
pub const __clockid_t_defined = @as(c_int, 1);
pub const __time_t_defined = @as(c_int, 1);
pub const __timer_t_defined = @as(c_int, 1);
pub const _BITS_STDINT_INTN_H = @as(c_int, 1);
pub const __BIT_TYPES_DEFINED__ = @as(c_int, 1);
pub const _ENDIAN_H = @as(c_int, 1);
pub const _BITS_ENDIAN_H = @as(c_int, 1);
pub const __LITTLE_ENDIAN = @as(c_int, 1234);
pub const __BIG_ENDIAN = @as(c_int, 4321);
pub const __PDP_ENDIAN = @as(c_int, 3412);
pub const _BITS_ENDIANNESS_H = @as(c_int, 1);
pub const __BYTE_ORDER = __LITTLE_ENDIAN;
pub const __FLOAT_WORD_ORDER = __BYTE_ORDER;
pub inline fn __LONG_LONG_PAIR(HI: anytype, LO: anytype) @TypeOf(HI) {
    _ = &HI;
    _ = &LO;
    return blk: {
        _ = &LO;
        break :blk HI;
    };
}
pub const LITTLE_ENDIAN = __LITTLE_ENDIAN;
pub const BIG_ENDIAN = __BIG_ENDIAN;
pub const PDP_ENDIAN = __PDP_ENDIAN;
pub const BYTE_ORDER = __BYTE_ORDER;
pub const _BITS_BYTESWAP_H = @as(c_int, 1);
pub inline fn __bswap_constant_16(x: anytype) __uint16_t {
    _ = &x;
    return @import("std").zig.c_translation.cast(__uint16_t, ((x >> @as(c_int, 8)) & @as(c_int, 0xff)) | ((x & @as(c_int, 0xff)) << @as(c_int, 8)));
}
pub inline fn __bswap_constant_32(x: anytype) @TypeOf(((((x & @import("std").zig.c_translation.promoteIntLiteral(c_uint, 0xff000000, .hex)) >> @as(c_int, 24)) | ((x & @import("std").zig.c_translation.promoteIntLiteral(c_uint, 0x00ff0000, .hex)) >> @as(c_int, 8))) | ((x & @as(c_uint, 0x0000ff00)) << @as(c_int, 8))) | ((x & @as(c_uint, 0x000000ff)) << @as(c_int, 24))) {
    _ = &x;
    return ((((x & @import("std").zig.c_translation.promoteIntLiteral(c_uint, 0xff000000, .hex)) >> @as(c_int, 24)) | ((x & @import("std").zig.c_translation.promoteIntLiteral(c_uint, 0x00ff0000, .hex)) >> @as(c_int, 8))) | ((x & @as(c_uint, 0x0000ff00)) << @as(c_int, 8))) | ((x & @as(c_uint, 0x000000ff)) << @as(c_int, 24));
}
pub inline fn __bswap_constant_64(x: anytype) @TypeOf(((((((((x & @as(c_ulonglong, 0xff00000000000000)) >> @as(c_int, 56)) | ((x & @as(c_ulonglong, 0x00ff000000000000)) >> @as(c_int, 40))) | ((x & @as(c_ulonglong, 0x0000ff0000000000)) >> @as(c_int, 24))) | ((x & @as(c_ulonglong, 0x000000ff00000000)) >> @as(c_int, 8))) | ((x & @as(c_ulonglong, 0x00000000ff000000)) << @as(c_int, 8))) | ((x & @as(c_ulonglong, 0x0000000000ff0000)) << @as(c_int, 24))) | ((x & @as(c_ulonglong, 0x000000000000ff00)) << @as(c_int, 40))) | ((x & @as(c_ulonglong, 0x00000000000000ff)) << @as(c_int, 56))) {
    _ = &x;
    return ((((((((x & @as(c_ulonglong, 0xff00000000000000)) >> @as(c_int, 56)) | ((x & @as(c_ulonglong, 0x00ff000000000000)) >> @as(c_int, 40))) | ((x & @as(c_ulonglong, 0x0000ff0000000000)) >> @as(c_int, 24))) | ((x & @as(c_ulonglong, 0x000000ff00000000)) >> @as(c_int, 8))) | ((x & @as(c_ulonglong, 0x00000000ff000000)) << @as(c_int, 8))) | ((x & @as(c_ulonglong, 0x0000000000ff0000)) << @as(c_int, 24))) | ((x & @as(c_ulonglong, 0x000000000000ff00)) << @as(c_int, 40))) | ((x & @as(c_ulonglong, 0x00000000000000ff)) << @as(c_int, 56));
}
pub const _BITS_UINTN_IDENTITY_H = @as(c_int, 1);
pub inline fn htobe16(x: anytype) @TypeOf(__bswap_16(x)) {
    _ = &x;
    return __bswap_16(x);
}
pub inline fn htole16(x: anytype) @TypeOf(__uint16_identity(x)) {
    _ = &x;
    return __uint16_identity(x);
}
pub inline fn be16toh(x: anytype) @TypeOf(__bswap_16(x)) {
    _ = &x;
    return __bswap_16(x);
}
pub inline fn le16toh(x: anytype) @TypeOf(__uint16_identity(x)) {
    _ = &x;
    return __uint16_identity(x);
}
pub inline fn htobe32(x: anytype) @TypeOf(__bswap_32(x)) {
    _ = &x;
    return __bswap_32(x);
}
pub inline fn htole32(x: anytype) @TypeOf(__uint32_identity(x)) {
    _ = &x;
    return __uint32_identity(x);
}
pub inline fn be32toh(x: anytype) @TypeOf(__bswap_32(x)) {
    _ = &x;
    return __bswap_32(x);
}
pub inline fn le32toh(x: anytype) @TypeOf(__uint32_identity(x)) {
    _ = &x;
    return __uint32_identity(x);
}
pub inline fn htobe64(x: anytype) @TypeOf(__bswap_64(x)) {
    _ = &x;
    return __bswap_64(x);
}
pub inline fn htole64(x: anytype) @TypeOf(__uint64_identity(x)) {
    _ = &x;
    return __uint64_identity(x);
}
pub inline fn be64toh(x: anytype) @TypeOf(__bswap_64(x)) {
    _ = &x;
    return __bswap_64(x);
}
pub inline fn le64toh(x: anytype) @TypeOf(__uint64_identity(x)) {
    _ = &x;
    return __uint64_identity(x);
}
pub const _SYS_SELECT_H = @as(c_int, 1);
pub const __FD_ZERO = @compileError("unable to translate macro: undefined identifier `__i`");
// /usr/include/bits/select.h:25:9
pub const __FD_SET = @compileError("unable to translate C expr: expected ')' instead got '|='");
// /usr/include/bits/select.h:32:9
pub const __FD_CLR = @compileError("unable to translate C expr: expected ')' instead got '&='");
// /usr/include/bits/select.h:34:9
pub inline fn __FD_ISSET(d: anytype, s: anytype) @TypeOf((__FDS_BITS(s)[@as(usize, @intCast(__FD_ELT(d)))] & __FD_MASK(d)) != @as(c_int, 0)) {
    _ = &d;
    _ = &s;
    return (__FDS_BITS(s)[@as(usize, @intCast(__FD_ELT(d)))] & __FD_MASK(d)) != @as(c_int, 0);
}
pub const __sigset_t_defined = @as(c_int, 1);
pub const ____sigset_t_defined = "";
pub const _SIGSET_NWORDS = @import("std").zig.c_translation.MacroArithmetic.div(@as(c_int, 1024), @as(c_int, 8) * @import("std").zig.c_translation.sizeof(c_ulong));
pub const __timeval_defined = @as(c_int, 1);
pub const _STRUCT_TIMESPEC = @as(c_int, 1);
pub const __suseconds_t_defined = "";
pub const __NFDBITS = @as(c_int, 8) * @import("std").zig.c_translation.cast(c_int, @import("std").zig.c_translation.sizeof(__fd_mask));
pub inline fn __FD_ELT(d: anytype) @TypeOf(@import("std").zig.c_translation.MacroArithmetic.div(d, __NFDBITS)) {
    _ = &d;
    return @import("std").zig.c_translation.MacroArithmetic.div(d, __NFDBITS);
}
pub inline fn __FD_MASK(d: anytype) __fd_mask {
    _ = &d;
    return @import("std").zig.c_translation.cast(__fd_mask, @as(c_ulong, 1) << @import("std").zig.c_translation.MacroArithmetic.rem(d, __NFDBITS));
}
pub inline fn __FDS_BITS(set: anytype) @TypeOf(set.*.__fds_bits) {
    _ = &set;
    return set.*.__fds_bits;
}
pub const FD_SETSIZE = __FD_SETSIZE;
pub const NFDBITS = __NFDBITS;
pub inline fn FD_SET(fd: anytype, fdsetp: anytype) @TypeOf(__FD_SET(fd, fdsetp)) {
    _ = &fd;
    _ = &fdsetp;
    return __FD_SET(fd, fdsetp);
}
pub inline fn FD_CLR(fd: anytype, fdsetp: anytype) @TypeOf(__FD_CLR(fd, fdsetp)) {
    _ = &fd;
    _ = &fdsetp;
    return __FD_CLR(fd, fdsetp);
}
pub inline fn FD_ISSET(fd: anytype, fdsetp: anytype) @TypeOf(__FD_ISSET(fd, fdsetp)) {
    _ = &fd;
    _ = &fdsetp;
    return __FD_ISSET(fd, fdsetp);
}
pub inline fn FD_ZERO(fdsetp: anytype) @TypeOf(__FD_ZERO(fdsetp)) {
    _ = &fdsetp;
    return __FD_ZERO(fdsetp);
}
pub const __blksize_t_defined = "";
pub const __blkcnt_t_defined = "";
pub const __fsblkcnt_t_defined = "";
pub const __fsfilcnt_t_defined = "";
pub const _BITS_PTHREADTYPES_COMMON_H = @as(c_int, 1);
pub const _THREAD_SHARED_TYPES_H = @as(c_int, 1);
pub const _BITS_PTHREADTYPES_ARCH_H = @as(c_int, 1);
pub const __SIZEOF_PTHREAD_MUTEX_T = @as(c_int, 40);
pub const __SIZEOF_PTHREAD_ATTR_T = @as(c_int, 56);
pub const __SIZEOF_PTHREAD_RWLOCK_T = @as(c_int, 56);
pub const __SIZEOF_PTHREAD_BARRIER_T = @as(c_int, 32);
pub const __SIZEOF_PTHREAD_MUTEXATTR_T = @as(c_int, 4);
pub const __SIZEOF_PTHREAD_COND_T = @as(c_int, 48);
pub const __SIZEOF_PTHREAD_CONDATTR_T = @as(c_int, 4);
pub const __SIZEOF_PTHREAD_RWLOCKATTR_T = @as(c_int, 8);
pub const __SIZEOF_PTHREAD_BARRIERATTR_T = @as(c_int, 4);
pub const __LOCK_ALIGNMENT = "";
pub const __ONCE_ALIGNMENT = "";
pub const _BITS_ATOMIC_WIDE_COUNTER_H = "";
pub const _THREAD_MUTEX_INTERNAL_H = @as(c_int, 1);
pub const __PTHREAD_MUTEX_HAVE_PREV = @as(c_int, 1);
pub const __PTHREAD_MUTEX_INITIALIZER = @compileError("unable to translate C expr: unexpected token '{'");
// /usr/include/bits/struct_mutex.h:56:10
pub const _RWLOCK_INTERNAL_H = "";
pub const __PTHREAD_RWLOCK_ELISION_EXTRA = @compileError("unable to translate C expr: unexpected token '{'");
// /usr/include/bits/struct_rwlock.h:40:11
pub inline fn __PTHREAD_RWLOCK_INITIALIZER(__flags: anytype) @TypeOf(__flags) {
    _ = &__flags;
    return blk: {
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = @as(c_int, 0);
        _ = &__PTHREAD_RWLOCK_ELISION_EXTRA;
        _ = @as(c_int, 0);
        break :blk __flags;
    };
}
pub const __ONCE_FLAG_INIT = @compileError("unable to translate C expr: unexpected token '{'");
// /usr/include/bits/thread-shared-types.h:112:9
pub const __have_pthread_attr_t = @as(c_int, 1);
pub const _ALLOCA_H = @as(c_int, 1);
pub const __COMPAR_FN_T = "";
pub const __GSL_MODE_H__ = "";
pub const __GSL_INLINE_H__ = "";
pub const INLINE_DECL = "";
pub inline fn GSL_RANGE_COND(x: anytype) @TypeOf(x) {
    _ = &x;
    return x;
}
pub const GSL_PREC_DOUBLE = @as(c_int, 0);
pub const GSL_PREC_SINGLE = @as(c_int, 1);
pub const GSL_PREC_APPROX = @as(c_int, 2);
pub inline fn GSL_MODE_PREC(mt: anytype) @TypeOf(mt & @import("std").zig.c_translation.cast(c_uint, @as(c_int, 7))) {
    _ = &mt;
    return mt & @import("std").zig.c_translation.cast(c_uint, @as(c_int, 7));
}
pub const GSL_MODE_DEFAULT = @as(c_int, 0);
pub const __GSL_PERMUTATION_H__ = "";
pub const __GSL_TYPES_H__ = "";
pub const GSL_VAR = @compileError("unable to translate C expr: unexpected token 'extern'");
// /usr/include/gsl/gsl_types.h:36:11
pub const __GSL_ERRNO_H__ = "";
pub const _STDIO_H = @as(c_int, 1);
pub const __need___va_list = "";
pub const __GNUC_VA_LIST = "";
pub const _____fpos_t_defined = @as(c_int, 1);
pub const ____mbstate_t_defined = @as(c_int, 1);
pub const _____fpos64_t_defined = @as(c_int, 1);
pub const ____FILE_defined = @as(c_int, 1);
pub const __FILE_defined = @as(c_int, 1);
pub const __struct_FILE_defined = @as(c_int, 1);
pub const __getc_unlocked_body = @compileError("TODO postfix inc/dec expr");
// /usr/include/bits/types/struct_FILE.h:105:9
pub const __putc_unlocked_body = @compileError("TODO postfix inc/dec expr");
// /usr/include/bits/types/struct_FILE.h:109:9
pub const _IO_EOF_SEEN = @as(c_int, 0x0010);
pub inline fn __feof_unlocked_body(_fp: anytype) @TypeOf((_fp.*._flags & _IO_EOF_SEEN) != @as(c_int, 0)) {
    _ = &_fp;
    return (_fp.*._flags & _IO_EOF_SEEN) != @as(c_int, 0);
}
pub const _IO_ERR_SEEN = @as(c_int, 0x0020);
pub inline fn __ferror_unlocked_body(_fp: anytype) @TypeOf((_fp.*._flags & _IO_ERR_SEEN) != @as(c_int, 0)) {
    _ = &_fp;
    return (_fp.*._flags & _IO_ERR_SEEN) != @as(c_int, 0);
}
pub const _IO_USER_LOCK = @import("std").zig.c_translation.promoteIntLiteral(c_int, 0x8000, .hex);
pub const __cookie_io_functions_t_defined = @as(c_int, 1);
pub const _VA_LIST_DEFINED = "";
pub const _IOFBF = @as(c_int, 0);
pub const _IOLBF = @as(c_int, 1);
pub const _IONBF = @as(c_int, 2);
pub const BUFSIZ = @as(c_int, 8192);
pub const EOF = -@as(c_int, 1);
pub const SEEK_SET = @as(c_int, 0);
pub const SEEK_CUR = @as(c_int, 1);
pub const SEEK_END = @as(c_int, 2);
pub const P_tmpdir = "/tmp";
pub const L_tmpnam = @as(c_int, 20);
pub const TMP_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 238328, .decimal);
pub const _BITS_STDIO_LIM_H = @as(c_int, 1);
pub const FILENAME_MAX = @as(c_int, 4096);
pub const L_ctermid = @as(c_int, 9);
pub const FOPEN_MAX = @as(c_int, 16);
pub const __attr_dealloc_fclose = __attr_dealloc(fclose, @as(c_int, 1));
pub const _ERRNO_H = @as(c_int, 1);
pub const _BITS_ERRNO_H = @as(c_int, 1);
pub const _ASM_GENERIC_ERRNO_H = "";
pub const _ASM_GENERIC_ERRNO_BASE_H = "";
pub const EPERM = @as(c_int, 1);
pub const ENOENT = @as(c_int, 2);
pub const ESRCH = @as(c_int, 3);
pub const EINTR = @as(c_int, 4);
pub const EIO = @as(c_int, 5);
pub const ENXIO = @as(c_int, 6);
pub const E2BIG = @as(c_int, 7);
pub const ENOEXEC = @as(c_int, 8);
pub const EBADF = @as(c_int, 9);
pub const ECHILD = @as(c_int, 10);
pub const EAGAIN = @as(c_int, 11);
pub const ENOMEM = @as(c_int, 12);
pub const EACCES = @as(c_int, 13);
pub const EFAULT = @as(c_int, 14);
pub const ENOTBLK = @as(c_int, 15);
pub const EBUSY = @as(c_int, 16);
pub const EEXIST = @as(c_int, 17);
pub const EXDEV = @as(c_int, 18);
pub const ENODEV = @as(c_int, 19);
pub const ENOTDIR = @as(c_int, 20);
pub const EISDIR = @as(c_int, 21);
pub const EINVAL = @as(c_int, 22);
pub const ENFILE = @as(c_int, 23);
pub const EMFILE = @as(c_int, 24);
pub const ENOTTY = @as(c_int, 25);
pub const ETXTBSY = @as(c_int, 26);
pub const EFBIG = @as(c_int, 27);
pub const ENOSPC = @as(c_int, 28);
pub const ESPIPE = @as(c_int, 29);
pub const EROFS = @as(c_int, 30);
pub const EMLINK = @as(c_int, 31);
pub const EPIPE = @as(c_int, 32);
pub const EDOM = @as(c_int, 33);
pub const ERANGE = @as(c_int, 34);
pub const EDEADLK = @as(c_int, 35);
pub const ENAMETOOLONG = @as(c_int, 36);
pub const ENOLCK = @as(c_int, 37);
pub const ENOSYS = @as(c_int, 38);
pub const ENOTEMPTY = @as(c_int, 39);
pub const ELOOP = @as(c_int, 40);
pub const EWOULDBLOCK = EAGAIN;
pub const ENOMSG = @as(c_int, 42);
pub const EIDRM = @as(c_int, 43);
pub const ECHRNG = @as(c_int, 44);
pub const EL2NSYNC = @as(c_int, 45);
pub const EL3HLT = @as(c_int, 46);
pub const EL3RST = @as(c_int, 47);
pub const ELNRNG = @as(c_int, 48);
pub const EUNATCH = @as(c_int, 49);
pub const ENOCSI = @as(c_int, 50);
pub const EL2HLT = @as(c_int, 51);
pub const EBADE = @as(c_int, 52);
pub const EBADR = @as(c_int, 53);
pub const EXFULL = @as(c_int, 54);
pub const ENOANO = @as(c_int, 55);
pub const EBADRQC = @as(c_int, 56);
pub const EBADSLT = @as(c_int, 57);
pub const EDEADLOCK = EDEADLK;
pub const EBFONT = @as(c_int, 59);
pub const ENOSTR = @as(c_int, 60);
pub const ENODATA = @as(c_int, 61);
pub const ETIME = @as(c_int, 62);
pub const ENOSR = @as(c_int, 63);
pub const ENONET = @as(c_int, 64);
pub const ENOPKG = @as(c_int, 65);
pub const EREMOTE = @as(c_int, 66);
pub const ENOLINK = @as(c_int, 67);
pub const EADV = @as(c_int, 68);
pub const ESRMNT = @as(c_int, 69);
pub const ECOMM = @as(c_int, 70);
pub const EPROTO = @as(c_int, 71);
pub const EMULTIHOP = @as(c_int, 72);
pub const EDOTDOT = @as(c_int, 73);
pub const EBADMSG = @as(c_int, 74);
pub const EOVERFLOW = @as(c_int, 75);
pub const ENOTUNIQ = @as(c_int, 76);
pub const EBADFD = @as(c_int, 77);
pub const EREMCHG = @as(c_int, 78);
pub const ELIBACC = @as(c_int, 79);
pub const ELIBBAD = @as(c_int, 80);
pub const ELIBSCN = @as(c_int, 81);
pub const ELIBMAX = @as(c_int, 82);
pub const ELIBEXEC = @as(c_int, 83);
pub const EILSEQ = @as(c_int, 84);
pub const ERESTART = @as(c_int, 85);
pub const ESTRPIPE = @as(c_int, 86);
pub const EUSERS = @as(c_int, 87);
pub const ENOTSOCK = @as(c_int, 88);
pub const EDESTADDRREQ = @as(c_int, 89);
pub const EMSGSIZE = @as(c_int, 90);
pub const EPROTOTYPE = @as(c_int, 91);
pub const ENOPROTOOPT = @as(c_int, 92);
pub const EPROTONOSUPPORT = @as(c_int, 93);
pub const ESOCKTNOSUPPORT = @as(c_int, 94);
pub const EOPNOTSUPP = @as(c_int, 95);
pub const EPFNOSUPPORT = @as(c_int, 96);
pub const EAFNOSUPPORT = @as(c_int, 97);
pub const EADDRINUSE = @as(c_int, 98);
pub const EADDRNOTAVAIL = @as(c_int, 99);
pub const ENETDOWN = @as(c_int, 100);
pub const ENETUNREACH = @as(c_int, 101);
pub const ENETRESET = @as(c_int, 102);
pub const ECONNABORTED = @as(c_int, 103);
pub const ECONNRESET = @as(c_int, 104);
pub const ENOBUFS = @as(c_int, 105);
pub const EISCONN = @as(c_int, 106);
pub const ENOTCONN = @as(c_int, 107);
pub const ESHUTDOWN = @as(c_int, 108);
pub const ETOOMANYREFS = @as(c_int, 109);
pub const ETIMEDOUT = @as(c_int, 110);
pub const ECONNREFUSED = @as(c_int, 111);
pub const EHOSTDOWN = @as(c_int, 112);
pub const EHOSTUNREACH = @as(c_int, 113);
pub const EALREADY = @as(c_int, 114);
pub const EINPROGRESS = @as(c_int, 115);
pub const ESTALE = @as(c_int, 116);
pub const EUCLEAN = @as(c_int, 117);
pub const ENOTNAM = @as(c_int, 118);
pub const ENAVAIL = @as(c_int, 119);
pub const EISNAM = @as(c_int, 120);
pub const EREMOTEIO = @as(c_int, 121);
pub const EDQUOT = @as(c_int, 122);
pub const ENOMEDIUM = @as(c_int, 123);
pub const EMEDIUMTYPE = @as(c_int, 124);
pub const ECANCELED = @as(c_int, 125);
pub const ENOKEY = @as(c_int, 126);
pub const EKEYEXPIRED = @as(c_int, 127);
pub const EKEYREVOKED = @as(c_int, 128);
pub const EKEYREJECTED = @as(c_int, 129);
pub const EOWNERDEAD = @as(c_int, 130);
pub const ENOTRECOVERABLE = @as(c_int, 131);
pub const ERFKILL = @as(c_int, 132);
pub const EHWPOISON = @as(c_int, 133);
pub const ENOTSUP = EOPNOTSUPP;
pub const errno = __errno_location().*;
pub const GSL_ERROR = @compileError("unable to translate macro: undefined identifier `__FILE__`");
// /usr/include/gsl/gsl_errno.h:104:9
pub const GSL_ERROR_VAL = @compileError("unable to translate macro: undefined identifier `__FILE__`");
// /usr/include/gsl/gsl_errno.h:112:9
pub const GSL_ERROR_VOID = @compileError("unable to translate macro: undefined identifier `__FILE__`");
// /usr/include/gsl/gsl_errno.h:121:9
pub inline fn GSL_ERROR_NULL(reason: anytype, gsl_errno: anytype) @TypeOf(GSL_ERROR_VAL(reason, gsl_errno, @as(c_int, 0))) {
    _ = &reason;
    _ = &gsl_errno;
    return GSL_ERROR_VAL(reason, gsl_errno, @as(c_int, 0));
}
pub inline fn GSL_ERROR_SELECT_2(a: anytype, b: anytype) @TypeOf(if (a != GSL_SUCCESS) a else if (b != GSL_SUCCESS) b else GSL_SUCCESS) {
    _ = &a;
    _ = &b;
    return if (a != GSL_SUCCESS) a else if (b != GSL_SUCCESS) b else GSL_SUCCESS;
}
pub inline fn GSL_ERROR_SELECT_3(a: anytype, b: anytype, c: anytype) @TypeOf(if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_2(b, c)) {
    _ = &a;
    _ = &b;
    _ = &c;
    return if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_2(b, c);
}
pub inline fn GSL_ERROR_SELECT_4(a: anytype, b: anytype, c: anytype, d: anytype) @TypeOf(if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_3(b, c, d)) {
    _ = &a;
    _ = &b;
    _ = &c;
    _ = &d;
    return if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_3(b, c, d);
}
pub inline fn GSL_ERROR_SELECT_5(a: anytype, b: anytype, c: anytype, d: anytype, e: anytype) @TypeOf(if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_4(b, c, d, e)) {
    _ = &a;
    _ = &b;
    _ = &c;
    _ = &d;
    _ = &e;
    return if (a != GSL_SUCCESS) a else GSL_ERROR_SELECT_4(b, c, d, e);
}
pub const GSL_STATUS_UPDATE = @compileError("unable to translate C expr: unexpected token 'do'");
// /usr/include/gsl/gsl_errno.h:150:9
pub const __GSL_CHECK_RANGE_H__ = "";
pub const GSL_RANGE_CHECK = @as(c_int, 1);
pub const __GSL_VECTOR_H__ = "";
pub const __GSL_VECTOR_COMPLEX_LONG_DOUBLE_H__ = "";
pub const __GSL_COMPLEX_H__ = "";
pub const GSL_COMPLEX_LEGACY = @as(c_int, 1);
pub const GSL_COMPLEX_DEFINE = @compileError("unable to translate macro: undefined identifier `dat`");
// /usr/include/gsl/gsl_complex.h:121:11
pub inline fn GSL_REAL(z: anytype) @TypeOf(z.dat[@as(usize, @intCast(@as(c_int, 0)))]) {
    _ = &z;
    return z.dat[@as(usize, @intCast(@as(c_int, 0)))];
}
pub inline fn GSL_IMAG(z: anytype) @TypeOf(z.dat[@as(usize, @intCast(@as(c_int, 1)))]) {
    _ = &z;
    return z.dat[@as(usize, @intCast(@as(c_int, 1)))];
}
pub inline fn GSL_COMPLEX_P(zp: anytype) @TypeOf(zp.*.dat) {
    _ = &zp;
    return zp.*.dat;
}
pub inline fn GSL_COMPLEX_P_REAL(zp: anytype) @TypeOf(zp.*.dat[@as(usize, @intCast(@as(c_int, 0)))]) {
    _ = &zp;
    return zp.*.dat[@as(usize, @intCast(@as(c_int, 0)))];
}
pub inline fn GSL_COMPLEX_P_IMAG(zp: anytype) @TypeOf(zp.*.dat[@as(usize, @intCast(@as(c_int, 1)))]) {
    _ = &zp;
    return zp.*.dat[@as(usize, @intCast(@as(c_int, 1)))];
}
pub inline fn GSL_COMPLEX_EQ(z1: anytype, z2: anytype) @TypeOf((z1.dat[@as(usize, @intCast(@as(c_int, 0)))] == z2.dat[@as(usize, @intCast(@as(c_int, 0)))]) and (z1.dat[@as(usize, @intCast(@as(c_int, 1)))] == z2.dat[@as(usize, @intCast(@as(c_int, 1)))])) {
    _ = &z1;
    _ = &z2;
    return (z1.dat[@as(usize, @intCast(@as(c_int, 0)))] == z2.dat[@as(usize, @intCast(@as(c_int, 0)))]) and (z1.dat[@as(usize, @intCast(@as(c_int, 1)))] == z2.dat[@as(usize, @intCast(@as(c_int, 1)))]);
}
pub const GSL_SET_COMPLEX = @compileError("unable to translate C expr: unexpected token 'do'");
// /usr/include/gsl/gsl_complex.h:130:11
pub const GSL_SET_REAL = @compileError("unable to translate C expr: unexpected token 'do'");
// /usr/include/gsl/gsl_complex.h:131:11
pub const GSL_SET_IMAG = @compileError("unable to translate C expr: unexpected token 'do'");
// /usr/include/gsl/gsl_complex.h:132:11
pub const GSL_SET_COMPLEX_PACKED = @compileError("unable to translate C expr: unexpected token 'do'");
// /usr/include/gsl/gsl_complex.h:140:9
pub const __GSL_VECTOR_LONG_DOUBLE_H__ = "";
pub const __GSL_BLOCK_LONG_DOUBLE_H__ = "";
pub const __GSL_VECTOR_COMPLEX_H__ = "";
pub inline fn GSL_VECTOR_REAL(z: anytype, i: anytype) @TypeOf(z.*.data[@as(usize, @intCast((@as(c_int, 2) * i) * z.*.stride))]) {
    _ = &z;
    _ = &i;
    return z.*.data[@as(usize, @intCast((@as(c_int, 2) * i) * z.*.stride))];
}
pub inline fn GSL_VECTOR_IMAG(z: anytype, i: anytype) @TypeOf(z.*.data[@as(usize, @intCast(((@as(c_int, 2) * i) * z.*.stride) + @as(c_int, 1)))]) {
    _ = &z;
    _ = &i;
    return z.*.data[@as(usize, @intCast(((@as(c_int, 2) * i) * z.*.stride) + @as(c_int, 1)))];
}
pub const GSL_VECTOR_COMPLEX = @compileError("unable to translate macro: undefined identifier `__FILE__`");
// /usr/include/gsl/gsl_vector_complex.h:8:9
pub inline fn GSL_COMPLEX_AT(zv: anytype, i: anytype) [*c]gsl_complex {
    _ = &zv;
    _ = &i;
    return @import("std").zig.c_translation.cast([*c]gsl_complex, &zv.*.data[@as(usize, @intCast((@as(c_int, 2) * i) * zv.*.stride))]);
}
pub inline fn GSL_COMPLEX_FLOAT_AT(zv: anytype, i: anytype) [*c]gsl_complex_float {
    _ = &zv;
    _ = &i;
    return @import("std").zig.c_translation.cast([*c]gsl_complex_float, &zv.*.data[@as(usize, @intCast((@as(c_int, 2) * i) * zv.*.stride))]);
}
pub inline fn GSL_COMPLEX_LONG_DOUBLE_AT(zv: anytype, i: anytype) [*c]gsl_complex_long_double {
    _ = &zv;
    _ = &i;
    return @import("std").zig.c_translation.cast([*c]gsl_complex_long_double, &zv.*.data[@as(usize, @intCast((@as(c_int, 2) * i) * zv.*.stride))]);
}
pub const __GSL_BLOCK_COMPLEX_LONG_DOUBLE_H__ = "";
pub const __GSL_VECTOR_COMPLEX_DOUBLE_H__ = "";
pub const __GSL_VECTOR_DOUBLE_H__ = "";
pub const __GSL_BLOCK_DOUBLE_H__ = "";
pub const __GSL_BLOCK_COMPLEX_DOUBLE_H__ = "";
pub const __GSL_VECTOR_COMPLEX_FLOAT_H__ = "";
pub const __GSL_VECTOR_FLOAT_H__ = "";
pub const __GSL_BLOCK_FLOAT_H__ = "";
pub const __GSL_BLOCK_COMPLEX_FLOAT_H__ = "";
pub const __GSL_VECTOR_ULONG_H__ = "";
pub const __GSL_BLOCK_ULONG_H__ = "";
pub const __GSL_VECTOR_LONG_H__ = "";
pub const __GSL_BLOCK_LONG_H__ = "";
pub const __GSL_VECTOR_UINT_H__ = "";
pub const __GSL_BLOCK_UINT_H__ = "";
pub const __GSL_VECTOR_INT_H__ = "";
pub const __GSL_BLOCK_INT_H__ = "";
pub const __GSL_VECTOR_USHORT_H__ = "";
pub const __GSL_BLOCK_USHORT_H__ = "";
pub const __GSL_VECTOR_SHORT_H__ = "";
pub const __GSL_BLOCK_SHORT_H__ = "";
pub const __GSL_VECTOR_UCHAR_H__ = "";
pub const __GSL_BLOCK_UCHAR_H__ = "";
pub const __GSL_VECTOR_CHAR_H__ = "";
pub const __GSL_BLOCK_CHAR_H__ = "";
pub const __GSL_MATRIX_H__ = "";
pub const __GSL_MATRIX_COMPLEX_LONG_DOUBLE_H__ = "";
pub const __GSL_BLAS_TYPES_H__ = "";
pub const __GSL_CBLAS_H__ = "";
pub const __need_ptrdiff_t = "";
pub const __need_max_align_t = "";
pub const __need_offsetof = "";
pub const __STDDEF_H = "";
pub const _PTRDIFF_T = "";
pub const __CLANG_MAX_ALIGN_T_DEFINED = "";
pub const offsetof = @compileError("unable to translate C expr: unexpected token 'an identifier'");
// /home/choros/zig/lib/include/__stddef_offsetof.h:16:9
pub const CBLAS_INDEX = usize;
pub const __GSL_MATRIX_COMPLEX_DOUBLE_H__ = "";
pub const __GSL_MATRIX_COMPLEX_FLOAT_H__ = "";
pub const __GSL_MATRIX_LONG_DOUBLE_H__ = "";
pub const __GSL_MATRIX_DOUBLE_H__ = "";
pub const __GSL_MATRIX_FLOAT_H__ = "";
pub const __GSL_MATRIX_ULONG_H__ = "";
pub const __GSL_MATRIX_LONG_H__ = "";
pub const __GSL_MATRIX_UINT_H__ = "";
pub const __GSL_MATRIX_INT_H__ = "";
pub const __GSL_MATRIX_USHORT_H__ = "";
pub const __GSL_MATRIX_SHORT_H__ = "";
pub const __GSL_MATRIX_UCHAR_H__ = "";
pub const __GSL_MATRIX_CHAR_H__ = "";
pub const __GSL_MATH_H__ = "";
pub const _MATH_H = @as(c_int, 1);
pub const _BITS_LIBM_SIMD_DECL_STUBS_H = @as(c_int, 1);
pub const __DECL_SIMD_cos = "";
pub const __DECL_SIMD_cosf = "";
pub const __DECL_SIMD_cosl = "";
pub const __DECL_SIMD_cosf16 = "";
pub const __DECL_SIMD_cosf32 = "";
pub const __DECL_SIMD_cosf64 = "";
pub const __DECL_SIMD_cosf128 = "";
pub const __DECL_SIMD_cosf32x = "";
pub const __DECL_SIMD_cosf64x = "";
pub const __DECL_SIMD_cosf128x = "";
pub const __DECL_SIMD_sin = "";
pub const __DECL_SIMD_sinf = "";
pub const __DECL_SIMD_sinl = "";
pub const __DECL_SIMD_sinf16 = "";
pub const __DECL_SIMD_sinf32 = "";
pub const __DECL_SIMD_sinf64 = "";
pub const __DECL_SIMD_sinf128 = "";
pub const __DECL_SIMD_sinf32x = "";
pub const __DECL_SIMD_sinf64x = "";
pub const __DECL_SIMD_sinf128x = "";
pub const __DECL_SIMD_sincos = "";
pub const __DECL_SIMD_sincosf = "";
pub const __DECL_SIMD_sincosl = "";
pub const __DECL_SIMD_sincosf16 = "";
pub const __DECL_SIMD_sincosf32 = "";
pub const __DECL_SIMD_sincosf64 = "";
pub const __DECL_SIMD_sincosf128 = "";
pub const __DECL_SIMD_sincosf32x = "";
pub const __DECL_SIMD_sincosf64x = "";
pub const __DECL_SIMD_sincosf128x = "";
pub const __DECL_SIMD_log = "";
pub const __DECL_SIMD_logf = "";
pub const __DECL_SIMD_logl = "";
pub const __DECL_SIMD_logf16 = "";
pub const __DECL_SIMD_logf32 = "";
pub const __DECL_SIMD_logf64 = "";
pub const __DECL_SIMD_logf128 = "";
pub const __DECL_SIMD_logf32x = "";
pub const __DECL_SIMD_logf64x = "";
pub const __DECL_SIMD_logf128x = "";
pub const __DECL_SIMD_exp = "";
pub const __DECL_SIMD_expf = "";
pub const __DECL_SIMD_expl = "";
pub const __DECL_SIMD_expf16 = "";
pub const __DECL_SIMD_expf32 = "";
pub const __DECL_SIMD_expf64 = "";
pub const __DECL_SIMD_expf128 = "";
pub const __DECL_SIMD_expf32x = "";
pub const __DECL_SIMD_expf64x = "";
pub const __DECL_SIMD_expf128x = "";
pub const __DECL_SIMD_pow = "";
pub const __DECL_SIMD_powf = "";
pub const __DECL_SIMD_powl = "";
pub const __DECL_SIMD_powf16 = "";
pub const __DECL_SIMD_powf32 = "";
pub const __DECL_SIMD_powf64 = "";
pub const __DECL_SIMD_powf128 = "";
pub const __DECL_SIMD_powf32x = "";
pub const __DECL_SIMD_powf64x = "";
pub const __DECL_SIMD_powf128x = "";
pub const __DECL_SIMD_acos = "";
pub const __DECL_SIMD_acosf = "";
pub const __DECL_SIMD_acosl = "";
pub const __DECL_SIMD_acosf16 = "";
pub const __DECL_SIMD_acosf32 = "";
pub const __DECL_SIMD_acosf64 = "";
pub const __DECL_SIMD_acosf128 = "";
pub const __DECL_SIMD_acosf32x = "";
pub const __DECL_SIMD_acosf64x = "";
pub const __DECL_SIMD_acosf128x = "";
pub const __DECL_SIMD_atan = "";
pub const __DECL_SIMD_atanf = "";
pub const __DECL_SIMD_atanl = "";
pub const __DECL_SIMD_atanf16 = "";
pub const __DECL_SIMD_atanf32 = "";
pub const __DECL_SIMD_atanf64 = "";
pub const __DECL_SIMD_atanf128 = "";
pub const __DECL_SIMD_atanf32x = "";
pub const __DECL_SIMD_atanf64x = "";
pub const __DECL_SIMD_atanf128x = "";
pub const __DECL_SIMD_asin = "";
pub const __DECL_SIMD_asinf = "";
pub const __DECL_SIMD_asinl = "";
pub const __DECL_SIMD_asinf16 = "";
pub const __DECL_SIMD_asinf32 = "";
pub const __DECL_SIMD_asinf64 = "";
pub const __DECL_SIMD_asinf128 = "";
pub const __DECL_SIMD_asinf32x = "";
pub const __DECL_SIMD_asinf64x = "";
pub const __DECL_SIMD_asinf128x = "";
pub const __DECL_SIMD_hypot = "";
pub const __DECL_SIMD_hypotf = "";
pub const __DECL_SIMD_hypotl = "";
pub const __DECL_SIMD_hypotf16 = "";
pub const __DECL_SIMD_hypotf32 = "";
pub const __DECL_SIMD_hypotf64 = "";
pub const __DECL_SIMD_hypotf128 = "";
pub const __DECL_SIMD_hypotf32x = "";
pub const __DECL_SIMD_hypotf64x = "";
pub const __DECL_SIMD_hypotf128x = "";
pub const __DECL_SIMD_exp2 = "";
pub const __DECL_SIMD_exp2f = "";
pub const __DECL_SIMD_exp2l = "";
pub const __DECL_SIMD_exp2f16 = "";
pub const __DECL_SIMD_exp2f32 = "";
pub const __DECL_SIMD_exp2f64 = "";
pub const __DECL_SIMD_exp2f128 = "";
pub const __DECL_SIMD_exp2f32x = "";
pub const __DECL_SIMD_exp2f64x = "";
pub const __DECL_SIMD_exp2f128x = "";
pub const __DECL_SIMD_exp10 = "";
pub const __DECL_SIMD_exp10f = "";
pub const __DECL_SIMD_exp10l = "";
pub const __DECL_SIMD_exp10f16 = "";
pub const __DECL_SIMD_exp10f32 = "";
pub const __DECL_SIMD_exp10f64 = "";
pub const __DECL_SIMD_exp10f128 = "";
pub const __DECL_SIMD_exp10f32x = "";
pub const __DECL_SIMD_exp10f64x = "";
pub const __DECL_SIMD_exp10f128x = "";
pub const __DECL_SIMD_cosh = "";
pub const __DECL_SIMD_coshf = "";
pub const __DECL_SIMD_coshl = "";
pub const __DECL_SIMD_coshf16 = "";
pub const __DECL_SIMD_coshf32 = "";
pub const __DECL_SIMD_coshf64 = "";
pub const __DECL_SIMD_coshf128 = "";
pub const __DECL_SIMD_coshf32x = "";
pub const __DECL_SIMD_coshf64x = "";
pub const __DECL_SIMD_coshf128x = "";
pub const __DECL_SIMD_expm1 = "";
pub const __DECL_SIMD_expm1f = "";
pub const __DECL_SIMD_expm1l = "";
pub const __DECL_SIMD_expm1f16 = "";
pub const __DECL_SIMD_expm1f32 = "";
pub const __DECL_SIMD_expm1f64 = "";
pub const __DECL_SIMD_expm1f128 = "";
pub const __DECL_SIMD_expm1f32x = "";
pub const __DECL_SIMD_expm1f64x = "";
pub const __DECL_SIMD_expm1f128x = "";
pub const __DECL_SIMD_sinh = "";
pub const __DECL_SIMD_sinhf = "";
pub const __DECL_SIMD_sinhl = "";
pub const __DECL_SIMD_sinhf16 = "";
pub const __DECL_SIMD_sinhf32 = "";
pub const __DECL_SIMD_sinhf64 = "";
pub const __DECL_SIMD_sinhf128 = "";
pub const __DECL_SIMD_sinhf32x = "";
pub const __DECL_SIMD_sinhf64x = "";
pub const __DECL_SIMD_sinhf128x = "";
pub const __DECL_SIMD_cbrt = "";
pub const __DECL_SIMD_cbrtf = "";
pub const __DECL_SIMD_cbrtl = "";
pub const __DECL_SIMD_cbrtf16 = "";
pub const __DECL_SIMD_cbrtf32 = "";
pub const __DECL_SIMD_cbrtf64 = "";
pub const __DECL_SIMD_cbrtf128 = "";
pub const __DECL_SIMD_cbrtf32x = "";
pub const __DECL_SIMD_cbrtf64x = "";
pub const __DECL_SIMD_cbrtf128x = "";
pub const __DECL_SIMD_atan2 = "";
pub const __DECL_SIMD_atan2f = "";
pub const __DECL_SIMD_atan2l = "";
pub const __DECL_SIMD_atan2f16 = "";
pub const __DECL_SIMD_atan2f32 = "";
pub const __DECL_SIMD_atan2f64 = "";
pub const __DECL_SIMD_atan2f128 = "";
pub const __DECL_SIMD_atan2f32x = "";
pub const __DECL_SIMD_atan2f64x = "";
pub const __DECL_SIMD_atan2f128x = "";
pub const __DECL_SIMD_log10 = "";
pub const __DECL_SIMD_log10f = "";
pub const __DECL_SIMD_log10l = "";
pub const __DECL_SIMD_log10f16 = "";
pub const __DECL_SIMD_log10f32 = "";
pub const __DECL_SIMD_log10f64 = "";
pub const __DECL_SIMD_log10f128 = "";
pub const __DECL_SIMD_log10f32x = "";
pub const __DECL_SIMD_log10f64x = "";
pub const __DECL_SIMD_log10f128x = "";
pub const __DECL_SIMD_log2 = "";
pub const __DECL_SIMD_log2f = "";
pub const __DECL_SIMD_log2l = "";
pub const __DECL_SIMD_log2f16 = "";
pub const __DECL_SIMD_log2f32 = "";
pub const __DECL_SIMD_log2f64 = "";
pub const __DECL_SIMD_log2f128 = "";
pub const __DECL_SIMD_log2f32x = "";
pub const __DECL_SIMD_log2f64x = "";
pub const __DECL_SIMD_log2f128x = "";
pub const __DECL_SIMD_log1p = "";
pub const __DECL_SIMD_log1pf = "";
pub const __DECL_SIMD_log1pl = "";
pub const __DECL_SIMD_log1pf16 = "";
pub const __DECL_SIMD_log1pf32 = "";
pub const __DECL_SIMD_log1pf64 = "";
pub const __DECL_SIMD_log1pf128 = "";
pub const __DECL_SIMD_log1pf32x = "";
pub const __DECL_SIMD_log1pf64x = "";
pub const __DECL_SIMD_log1pf128x = "";
pub const __DECL_SIMD_logp1 = "";
pub const __DECL_SIMD_logp1f = "";
pub const __DECL_SIMD_logp1l = "";
pub const __DECL_SIMD_logp1f16 = "";
pub const __DECL_SIMD_logp1f32 = "";
pub const __DECL_SIMD_logp1f64 = "";
pub const __DECL_SIMD_logp1f128 = "";
pub const __DECL_SIMD_logp1f32x = "";
pub const __DECL_SIMD_logp1f64x = "";
pub const __DECL_SIMD_logp1f128x = "";
pub const __DECL_SIMD_atanh = "";
pub const __DECL_SIMD_atanhf = "";
pub const __DECL_SIMD_atanhl = "";
pub const __DECL_SIMD_atanhf16 = "";
pub const __DECL_SIMD_atanhf32 = "";
pub const __DECL_SIMD_atanhf64 = "";
pub const __DECL_SIMD_atanhf128 = "";
pub const __DECL_SIMD_atanhf32x = "";
pub const __DECL_SIMD_atanhf64x = "";
pub const __DECL_SIMD_atanhf128x = "";
pub const __DECL_SIMD_acosh = "";
pub const __DECL_SIMD_acoshf = "";
pub const __DECL_SIMD_acoshl = "";
pub const __DECL_SIMD_acoshf16 = "";
pub const __DECL_SIMD_acoshf32 = "";
pub const __DECL_SIMD_acoshf64 = "";
pub const __DECL_SIMD_acoshf128 = "";
pub const __DECL_SIMD_acoshf32x = "";
pub const __DECL_SIMD_acoshf64x = "";
pub const __DECL_SIMD_acoshf128x = "";
pub const __DECL_SIMD_erf = "";
pub const __DECL_SIMD_erff = "";
pub const __DECL_SIMD_erfl = "";
pub const __DECL_SIMD_erff16 = "";
pub const __DECL_SIMD_erff32 = "";
pub const __DECL_SIMD_erff64 = "";
pub const __DECL_SIMD_erff128 = "";
pub const __DECL_SIMD_erff32x = "";
pub const __DECL_SIMD_erff64x = "";
pub const __DECL_SIMD_erff128x = "";
pub const __DECL_SIMD_tanh = "";
pub const __DECL_SIMD_tanhf = "";
pub const __DECL_SIMD_tanhl = "";
pub const __DECL_SIMD_tanhf16 = "";
pub const __DECL_SIMD_tanhf32 = "";
pub const __DECL_SIMD_tanhf64 = "";
pub const __DECL_SIMD_tanhf128 = "";
pub const __DECL_SIMD_tanhf32x = "";
pub const __DECL_SIMD_tanhf64x = "";
pub const __DECL_SIMD_tanhf128x = "";
pub const __DECL_SIMD_asinh = "";
pub const __DECL_SIMD_asinhf = "";
pub const __DECL_SIMD_asinhl = "";
pub const __DECL_SIMD_asinhf16 = "";
pub const __DECL_SIMD_asinhf32 = "";
pub const __DECL_SIMD_asinhf64 = "";
pub const __DECL_SIMD_asinhf128 = "";
pub const __DECL_SIMD_asinhf32x = "";
pub const __DECL_SIMD_asinhf64x = "";
pub const __DECL_SIMD_asinhf128x = "";
pub const __DECL_SIMD_erfc = "";
pub const __DECL_SIMD_erfcf = "";
pub const __DECL_SIMD_erfcl = "";
pub const __DECL_SIMD_erfcf16 = "";
pub const __DECL_SIMD_erfcf32 = "";
pub const __DECL_SIMD_erfcf64 = "";
pub const __DECL_SIMD_erfcf128 = "";
pub const __DECL_SIMD_erfcf32x = "";
pub const __DECL_SIMD_erfcf64x = "";
pub const __DECL_SIMD_erfcf128x = "";
pub const __DECL_SIMD_tan = "";
pub const __DECL_SIMD_tanf = "";
pub const __DECL_SIMD_tanl = "";
pub const __DECL_SIMD_tanf16 = "";
pub const __DECL_SIMD_tanf32 = "";
pub const __DECL_SIMD_tanf64 = "";
pub const __DECL_SIMD_tanf128 = "";
pub const __DECL_SIMD_tanf32x = "";
pub const __DECL_SIMD_tanf64x = "";
pub const __DECL_SIMD_tanf128x = "";
pub const __DECL_SIMD_sinpi = "";
pub const __DECL_SIMD_sinpif = "";
pub const __DECL_SIMD_sinpil = "";
pub const __DECL_SIMD_sinpif16 = "";
pub const __DECL_SIMD_sinpif32 = "";
pub const __DECL_SIMD_sinpif64 = "";
pub const __DECL_SIMD_sinpif128 = "";
pub const __DECL_SIMD_sinpif32x = "";
pub const __DECL_SIMD_sinpif64x = "";
pub const __DECL_SIMD_sinpif128x = "";
pub const __DECL_SIMD_cospi = "";
pub const __DECL_SIMD_cospif = "";
pub const __DECL_SIMD_cospil = "";
pub const __DECL_SIMD_cospif16 = "";
pub const __DECL_SIMD_cospif32 = "";
pub const __DECL_SIMD_cospif64 = "";
pub const __DECL_SIMD_cospif128 = "";
pub const __DECL_SIMD_cospif32x = "";
pub const __DECL_SIMD_cospif64x = "";
pub const __DECL_SIMD_cospif128x = "";
pub const __DECL_SIMD_tanpi = "";
pub const __DECL_SIMD_tanpif = "";
pub const __DECL_SIMD_tanpil = "";
pub const __DECL_SIMD_tanpif16 = "";
pub const __DECL_SIMD_tanpif32 = "";
pub const __DECL_SIMD_tanpif64 = "";
pub const __DECL_SIMD_tanpif128 = "";
pub const __DECL_SIMD_tanpif32x = "";
pub const __DECL_SIMD_tanpif64x = "";
pub const __DECL_SIMD_tanpif128x = "";
pub const HUGE_VAL = @compileError("unable to translate macro: undefined identifier `__builtin_huge_val`");
// /usr/include/math.h:48:10
pub const HUGE_VALF = __builtin_huge_valf();
pub const HUGE_VALL = @compileError("unable to translate macro: undefined identifier `__builtin_huge_vall`");
// /usr/include/math.h:60:11
pub const INFINITY = __builtin_inff();
pub const NAN = __builtin_nanf("");
pub const __GLIBC_FLT_EVAL_METHOD = @compileError("unable to translate macro: undefined identifier `__FLT_EVAL_METHOD__`");
// /usr/include/bits/flt-eval-method.h:27:11
pub const __FP_LOGB0_IS_MIN = @as(c_int, 1);
pub const __FP_LOGBNAN_IS_MIN = @as(c_int, 1);
pub const FP_ILOGB0 = -@import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal) - @as(c_int, 1);
pub const FP_ILOGBNAN = -@import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal) - @as(c_int, 1);
pub const __SIMD_DECL = @compileError("unable to translate macro: undefined identifier `__DECL_SIMD_`");
// /usr/include/bits/mathcalls-macros.h:19:9
pub inline fn __MATHCALL_VEC(function: anytype, suffix: anytype, args: anytype) @TypeOf(__SIMD_DECL(__MATH_PRECNAME(function, suffix)) ++ __MATHCALL(function, suffix, args)) {
    _ = &function;
    _ = &suffix;
    _ = &args;
    return __SIMD_DECL(__MATH_PRECNAME(function, suffix)) ++ __MATHCALL(function, suffix, args);
}
pub inline fn __MATHDECL_VEC(@"type": anytype, function: anytype, suffix: anytype, args: anytype) @TypeOf(__SIMD_DECL(__MATH_PRECNAME(function, suffix)) ++ __MATHDECL(@"type", function, suffix, args)) {
    _ = &@"type";
    _ = &function;
    _ = &suffix;
    _ = &args;
    return __SIMD_DECL(__MATH_PRECNAME(function, suffix)) ++ __MATHDECL(@"type", function, suffix, args);
}
pub inline fn __MATHCALL(function: anytype, suffix: anytype, args: anytype) @TypeOf(__MATHDECL(_Mdouble_, function, suffix, args)) {
    _ = &function;
    _ = &suffix;
    _ = &args;
    return __MATHDECL(_Mdouble_, function, suffix, args);
}
pub const __MATHDECL = @compileError("unable to translate macro: undefined identifier `__`");
// /usr/include/bits/mathcalls-macros.h:31:9
pub inline fn __MATHCALLX(function: anytype, suffix: anytype, args: anytype, attrib: anytype) @TypeOf(__MATHDECLX(_Mdouble_, function, suffix, args, attrib)) {
    _ = &function;
    _ = &suffix;
    _ = &args;
    _ = &attrib;
    return __MATHDECLX(_Mdouble_, function, suffix, args, attrib);
}
pub const __MATHDECLX = @compileError("unable to translate C expr: unexpected token '__attribute__'");
// /usr/include/bits/mathcalls-macros.h:36:9
pub const __MATHDECL_1_IMPL = @compileError("unable to translate C expr: unexpected token 'extern'");
// /usr/include/bits/mathcalls-macros.h:38:9
pub inline fn __MATHDECL_1(@"type": anytype, function: anytype, suffix: anytype, args: anytype) @TypeOf(__MATHDECL_1_IMPL(@"type", function, suffix, args)) {
    _ = &@"type";
    _ = &function;
    _ = &suffix;
    _ = &args;
    return __MATHDECL_1_IMPL(@"type", function, suffix, args);
}
pub inline fn __MATHDECL_ALIAS(@"type": anytype, function: anytype, suffix: anytype, args: anytype, alias: anytype) @TypeOf(__MATHDECL_1(@"type", function, suffix, args)) {
    _ = &@"type";
    _ = &function;
    _ = &suffix;
    _ = &args;
    _ = &alias;
    return __MATHDECL_1(@"type", function, suffix, args);
}
pub const __MATHREDIR = @compileError("unable to translate C expr: unexpected token 'extern'");
// /usr/include/bits/mathcalls-macros.h:47:9
pub const _Mdouble_ = f64;
pub inline fn __MATH_PRECNAME(name: anytype, r: anytype) @TypeOf(__CONCAT(name, r)) {
    _ = &name;
    _ = &r;
    return __CONCAT(name, r);
}
pub const __MATH_DECLARING_DOUBLE = @as(c_int, 1);
pub const __MATH_DECLARING_FLOATN = @as(c_int, 0);
pub const __MATH_DECLARE_LDOUBLE = @as(c_int, 1);
pub const __MATHCALL_NARROW_ARGS_1 = @compileError("unable to translate macro: undefined identifier `_Marg_`");
// /usr/include/math.h:519:9
pub const __MATHCALL_NARROW_ARGS_2 = @compileError("unable to translate macro: undefined identifier `_Marg_`");
// /usr/include/math.h:520:9
pub const __MATHCALL_NARROW_ARGS_3 = @compileError("unable to translate macro: undefined identifier `_Marg_`");
// /usr/include/math.h:521:9
pub const __MATHCALL_NARROW_NORMAL = @compileError("unable to translate macro: undefined identifier `_Mret_`");
// /usr/include/math.h:522:9
pub const __MATHCALL_NARROW_REDIR = @compileError("unable to translate macro: undefined identifier `_Mret_`");
// /usr/include/math.h:524:9
pub inline fn __MATHCALL_NARROW(func: anytype, redir: anytype, nargs: anytype) @TypeOf(__MATHCALL_NARROW_NORMAL(func, nargs)) {
    _ = &func;
    _ = &redir;
    _ = &nargs;
    return __MATHCALL_NARROW_NORMAL(func, nargs);
}
pub const __MATH_TG_F32 = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/math.h:855:12
pub const __MATH_TG_F64X = @compileError("unable to translate C expr: unexpected token ''");
// /usr/include/math.h:864:12
pub const __MATH_TG = @compileError("unable to translate macro: undefined identifier `f`");
// /usr/include/math.h:866:11
pub const fpclassify = @compileError("unable to translate macro: undefined identifier `__builtin_fpclassify`");
// /usr/include/math.h:936:11
pub inline fn signbit(x: anytype) @TypeOf(__builtin_signbit(x)) {
    _ = &x;
    return __builtin_signbit(x);
}
pub const isfinite = @compileError("unable to translate macro: undefined identifier `__builtin_isfinite`");
// /usr/include/math.h:963:11
pub const isnormal = @compileError("unable to translate macro: undefined identifier `__builtin_isnormal`");
// /usr/include/math.h:971:11
pub const MATH_ERRNO = @as(c_int, 1);
pub const MATH_ERREXCEPT = @as(c_int, 2);
pub const math_errhandling = MATH_ERRNO | MATH_ERREXCEPT;
pub const M_E = @as(f64, 2.7182818284590452354);
pub const M_LOG2E = @as(f64, 1.4426950408889634074);
pub const M_LOG10E = @as(f64, 0.43429448190325182765);
pub const M_LN2 = @as(f64, 0.69314718055994530942);
pub const M_LN10 = @as(f64, 2.30258509299404568402);
pub const M_PI = @as(f64, 3.14159265358979323846);
pub const M_PI_2 = @as(f64, 1.57079632679489661923);
pub const M_PI_4 = @as(f64, 0.78539816339744830962);
pub const M_1_PI = @as(f64, 0.31830988618379067154);
pub const M_2_PI = @as(f64, 0.63661977236758134308);
pub const M_2_SQRTPI = @as(f64, 1.12837916709551257390);
pub const M_SQRT2 = @as(f64, 1.41421356237309504880);
pub const M_SQRT1_2 = @as(f64, 0.70710678118654752440);
pub const isgreater = @compileError("unable to translate macro: undefined identifier `__builtin_isgreater`");
// /usr/include/math.h:1275:11
pub const isgreaterequal = @compileError("unable to translate macro: undefined identifier `__builtin_isgreaterequal`");
// /usr/include/math.h:1276:11
pub const isless = @compileError("unable to translate macro: undefined identifier `__builtin_isless`");
// /usr/include/math.h:1277:11
pub const islessequal = @compileError("unable to translate macro: undefined identifier `__builtin_islessequal`");
// /usr/include/math.h:1278:11
pub const islessgreater = @compileError("unable to translate macro: undefined identifier `__builtin_islessgreater`");
// /usr/include/math.h:1279:11
pub const isunordered = @compileError("unable to translate macro: undefined identifier `__builtin_isunordered`");
// /usr/include/math.h:1280:11
pub const __GSL_SYS_H__ = "";
pub const __GSL_MACHINE_H__ = "";
pub const __CLANG_LIMITS_H = "";
pub const _GCC_LIMITS_H_ = "";
pub const _LIBC_LIMITS_H_ = @as(c_int, 1);
pub const MB_LEN_MAX = @as(c_int, 16);
pub const LLONG_MIN = -LLONG_MAX - @as(c_int, 1);
pub const LLONG_MAX = __LONG_LONG_MAX__;
pub const ULLONG_MAX = (LLONG_MAX * @as(c_ulonglong, 2)) + @as(c_int, 1);
pub const _BITS_POSIX1_LIM_H = @as(c_int, 1);
pub const _POSIX_AIO_LISTIO_MAX = @as(c_int, 2);
pub const _POSIX_AIO_MAX = @as(c_int, 1);
pub const _POSIX_ARG_MAX = @as(c_int, 4096);
pub const _POSIX_CHILD_MAX = @as(c_int, 25);
pub const _POSIX_DELAYTIMER_MAX = @as(c_int, 32);
pub const _POSIX_HOST_NAME_MAX = @as(c_int, 255);
pub const _POSIX_LINK_MAX = @as(c_int, 8);
pub const _POSIX_LOGIN_NAME_MAX = @as(c_int, 9);
pub const _POSIX_MAX_CANON = @as(c_int, 255);
pub const _POSIX_MAX_INPUT = @as(c_int, 255);
pub const _POSIX_MQ_OPEN_MAX = @as(c_int, 8);
pub const _POSIX_MQ_PRIO_MAX = @as(c_int, 32);
pub const _POSIX_NAME_MAX = @as(c_int, 14);
pub const _POSIX_NGROUPS_MAX = @as(c_int, 8);
pub const _POSIX_OPEN_MAX = @as(c_int, 20);
pub const _POSIX_PATH_MAX = @as(c_int, 256);
pub const _POSIX_PIPE_BUF = @as(c_int, 512);
pub const _POSIX_RE_DUP_MAX = @as(c_int, 255);
pub const _POSIX_RTSIG_MAX = @as(c_int, 8);
pub const _POSIX_SEM_NSEMS_MAX = @as(c_int, 256);
pub const _POSIX_SEM_VALUE_MAX = @as(c_int, 32767);
pub const _POSIX_SIGQUEUE_MAX = @as(c_int, 32);
pub const _POSIX_SSIZE_MAX = @as(c_int, 32767);
pub const _POSIX_STREAM_MAX = @as(c_int, 8);
pub const _POSIX_SYMLINK_MAX = @as(c_int, 255);
pub const _POSIX_SYMLOOP_MAX = @as(c_int, 8);
pub const _POSIX_TIMER_MAX = @as(c_int, 32);
pub const _POSIX_TTY_NAME_MAX = @as(c_int, 9);
pub const _POSIX_TZNAME_MAX = @as(c_int, 6);
pub const _POSIX_CLOCKRES_MIN = @import("std").zig.c_translation.promoteIntLiteral(c_int, 20000000, .decimal);
pub const __undef_NR_OPEN = "";
pub const __undef_LINK_MAX = "";
pub const __undef_OPEN_MAX = "";
pub const __undef_ARG_MAX = "";
pub const _LINUX_LIMITS_H = "";
pub const NR_OPEN = @as(c_int, 1024);
pub const NGROUPS_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65536, .decimal);
pub const ARG_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 131072, .decimal);
pub const LINK_MAX = @as(c_int, 127);
pub const MAX_CANON = @as(c_int, 255);
pub const MAX_INPUT = @as(c_int, 255);
pub const NAME_MAX = @as(c_int, 255);
pub const PATH_MAX = @as(c_int, 4096);
pub const PIPE_BUF = @as(c_int, 4096);
pub const XATTR_NAME_MAX = @as(c_int, 255);
pub const XATTR_SIZE_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65536, .decimal);
pub const XATTR_LIST_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 65536, .decimal);
pub const RTSIG_MAX = @as(c_int, 32);
pub const _POSIX_THREAD_KEYS_MAX = @as(c_int, 128);
pub const PTHREAD_KEYS_MAX = @as(c_int, 1024);
pub const _POSIX_THREAD_DESTRUCTOR_ITERATIONS = @as(c_int, 4);
pub const PTHREAD_DESTRUCTOR_ITERATIONS = _POSIX_THREAD_DESTRUCTOR_ITERATIONS;
pub const _POSIX_THREAD_THREADS_MAX = @as(c_int, 64);
pub const AIO_PRIO_DELTA_MAX = @as(c_int, 20);
pub const PTHREAD_STACK_MIN = @as(c_int, 16384);
pub const DELAYTIMER_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const TTY_NAME_MAX = @as(c_int, 32);
pub const LOGIN_NAME_MAX = @as(c_int, 256);
pub const HOST_NAME_MAX = @as(c_int, 64);
pub const MQ_PRIO_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 32768, .decimal);
pub const SEM_VALUE_MAX = @import("std").zig.c_translation.promoteIntLiteral(c_int, 2147483647, .decimal);
pub const SSIZE_MAX = LONG_MAX;
pub const _BITS_POSIX2_LIM_H = @as(c_int, 1);
pub const _POSIX2_BC_BASE_MAX = @as(c_int, 99);
pub const _POSIX2_BC_DIM_MAX = @as(c_int, 2048);
pub const _POSIX2_BC_SCALE_MAX = @as(c_int, 99);
pub const _POSIX2_BC_STRING_MAX = @as(c_int, 1000);
pub const _POSIX2_COLL_WEIGHTS_MAX = @as(c_int, 2);
pub const _POSIX2_EXPR_NEST_MAX = @as(c_int, 32);
pub const _POSIX2_LINE_MAX = @as(c_int, 2048);
pub const _POSIX2_RE_DUP_MAX = @as(c_int, 255);
pub const _POSIX2_CHARCLASS_NAME_MAX = @as(c_int, 14);
pub const BC_BASE_MAX = _POSIX2_BC_BASE_MAX;
pub const BC_DIM_MAX = _POSIX2_BC_DIM_MAX;
pub const BC_SCALE_MAX = _POSIX2_BC_SCALE_MAX;
pub const BC_STRING_MAX = _POSIX2_BC_STRING_MAX;
pub const COLL_WEIGHTS_MAX = @as(c_int, 255);
pub const EXPR_NEST_MAX = _POSIX2_EXPR_NEST_MAX;
pub const LINE_MAX = _POSIX2_LINE_MAX;
pub const CHARCLASS_NAME_MAX = @as(c_int, 2048);
pub const RE_DUP_MAX = @as(c_int, 0x7fff);
pub const SCHAR_MAX = __SCHAR_MAX__;
pub const SHRT_MAX = __SHRT_MAX__;
pub const INT_MAX = __INT_MAX__;
pub const LONG_MAX = __LONG_MAX__;
pub const SCHAR_MIN = -__SCHAR_MAX__ - @as(c_int, 1);
pub const SHRT_MIN = -__SHRT_MAX__ - @as(c_int, 1);
pub const INT_MIN = -__INT_MAX__ - @as(c_int, 1);
pub const LONG_MIN = -__LONG_MAX__ - @as(c_long, 1);
pub const UCHAR_MAX = (__SCHAR_MAX__ * @as(c_int, 2)) + @as(c_int, 1);
pub const USHRT_MAX = (__SHRT_MAX__ * @as(c_int, 2)) + @as(c_int, 1);
pub const UINT_MAX = (__INT_MAX__ * @as(c_uint, 2)) + @as(c_uint, 1);
pub const ULONG_MAX = (__LONG_MAX__ * @as(c_ulong, 2)) + @as(c_ulong, 1);
pub const CHAR_BIT = __CHAR_BIT__;
pub const CHAR_MIN = SCHAR_MIN;
pub const CHAR_MAX = __SCHAR_MAX__;
pub const __CLANG_FLOAT_H = "";
pub const FLT_EVAL_METHOD = @compileError("unable to translate macro: undefined identifier `__FLT_EVAL_METHOD__`");
// /home/choros/zig/lib/include/float.h:107:9
pub const FLT_ROUNDS = @compileError("unable to translate macro: undefined identifier `__builtin_flt_rounds`");
// /home/choros/zig/lib/include/float.h:109:9
pub const FLT_RADIX = __FLT_RADIX__;
pub const FLT_MANT_DIG = __FLT_MANT_DIG__;
pub const DBL_MANT_DIG = __DBL_MANT_DIG__;
pub const LDBL_MANT_DIG = __LDBL_MANT_DIG__;
pub const DECIMAL_DIG = __DECIMAL_DIG__;
pub const FLT_DIG = __FLT_DIG__;
pub const DBL_DIG = __DBL_DIG__;
pub const LDBL_DIG = __LDBL_DIG__;
pub const FLT_MIN_EXP = __FLT_MIN_EXP__;
pub const DBL_MIN_EXP = __DBL_MIN_EXP__;
pub const LDBL_MIN_EXP = __LDBL_MIN_EXP__;
pub const FLT_MIN_10_EXP = __FLT_MIN_10_EXP__;
pub const DBL_MIN_10_EXP = __DBL_MIN_10_EXP__;
pub const LDBL_MIN_10_EXP = __LDBL_MIN_10_EXP__;
pub const FLT_MAX_EXP = __FLT_MAX_EXP__;
pub const DBL_MAX_EXP = __DBL_MAX_EXP__;
pub const LDBL_MAX_EXP = __LDBL_MAX_EXP__;
pub const FLT_MAX_10_EXP = __FLT_MAX_10_EXP__;
pub const DBL_MAX_10_EXP = __DBL_MAX_10_EXP__;
pub const LDBL_MAX_10_EXP = __LDBL_MAX_10_EXP__;
pub const FLT_MAX = __FLT_MAX__;
pub const DBL_MAX = __DBL_MAX__;
pub const LDBL_MAX = __LDBL_MAX__;
pub const FLT_EPSILON = __FLT_EPSILON__;
pub const DBL_EPSILON = __DBL_EPSILON__;
pub const LDBL_EPSILON = __LDBL_EPSILON__;
pub const FLT_MIN = __FLT_MIN__;
pub const DBL_MIN = __DBL_MIN__;
pub const LDBL_MIN = __LDBL_MIN__;
pub const FLT_TRUE_MIN = __FLT_DENORM_MIN__;
pub const DBL_TRUE_MIN = __DBL_DENORM_MIN__;
pub const LDBL_TRUE_MIN = __LDBL_DENORM_MIN__;
pub const FLT_DECIMAL_DIG = __FLT_DECIMAL_DIG__;
pub const DBL_DECIMAL_DIG = __DBL_DECIMAL_DIG__;
pub const LDBL_DECIMAL_DIG = __LDBL_DECIMAL_DIG__;
pub const FLT_HAS_SUBNORM = __FLT_HAS_DENORM__;
pub const DBL_HAS_SUBNORM = __DBL_HAS_DENORM__;
pub const LDBL_HAS_SUBNORM = __LDBL_HAS_DENORM__;
pub const FLT_NORM_MAX = __FLT_NORM_MAX__;
pub const DBL_NORM_MAX = __DBL_NORM_MAX__;
pub const LDBL_NORM_MAX = __LDBL_NORM_MAX__;
pub const GSL_DBL_EPSILON = @as(f64, 2.2204460492503131e-16);
pub const GSL_SQRT_DBL_EPSILON = @as(f64, 1.4901161193847656e-08);
pub const GSL_ROOT3_DBL_EPSILON = @as(f64, 6.0554544523933429e-06);
pub const GSL_ROOT4_DBL_EPSILON = @as(f64, 1.2207031250000000e-04);
pub const GSL_ROOT5_DBL_EPSILON = @as(f64, 7.4009597974140505e-04);
pub const GSL_ROOT6_DBL_EPSILON = @as(f64, 2.4607833005759251e-03);
pub const GSL_LOG_DBL_EPSILON = -@as(f64, 3.6043653389117154e+01);
pub const GSL_DBL_MIN = @as(f64, 2.2250738585072014e-308);
pub const GSL_SQRT_DBL_MIN = @as(f64, 1.4916681462400413e-154);
pub const GSL_ROOT3_DBL_MIN = @as(f64, 2.8126442852362996e-103);
pub const GSL_ROOT4_DBL_MIN = @as(f64, 1.2213386697554620e-77);
pub const GSL_ROOT5_DBL_MIN = @as(f64, 2.9476022969691763e-62);
pub const GSL_ROOT6_DBL_MIN = @as(f64, 5.3034368905798218e-52);
pub const GSL_LOG_DBL_MIN = -@as(f64, 7.0839641853226408e+02);
pub const GSL_DBL_MAX = @as(f64, 1.7976931348623157e+308);
pub const GSL_SQRT_DBL_MAX = @as(f64, 1.3407807929942596e+154);
pub const GSL_ROOT3_DBL_MAX = @as(f64, 5.6438030941222897e+102);
pub const GSL_ROOT4_DBL_MAX = @as(f64, 1.1579208923731620e+77);
pub const GSL_ROOT5_DBL_MAX = @as(f64, 4.4765466227572707e+61);
pub const GSL_ROOT6_DBL_MAX = @as(f64, 2.3756689782295612e+51);
pub const GSL_LOG_DBL_MAX = @as(f64, 7.0978271289338397e+02);
pub const GSL_FLT_EPSILON = @as(f64, 1.1920928955078125e-07);
pub const GSL_SQRT_FLT_EPSILON = @as(f64, 3.4526698300124393e-04);
pub const GSL_ROOT3_FLT_EPSILON = @as(f64, 4.9215666011518501e-03);
pub const GSL_ROOT4_FLT_EPSILON = @as(f64, 1.8581361171917516e-02);
pub const GSL_ROOT5_FLT_EPSILON = @as(f64, 4.1234622211652937e-02);
pub const GSL_ROOT6_FLT_EPSILON = @as(f64, 7.0153878019335827e-02);
pub const GSL_LOG_FLT_EPSILON = -@as(f64, 1.5942385152878742e+01);
pub const GSL_FLT_MIN = @as(f64, 1.1754943508222875e-38);
pub const GSL_SQRT_FLT_MIN = @as(f64, 1.0842021724855044e-19);
pub const GSL_ROOT3_FLT_MIN = @as(f64, 2.2737367544323241e-13);
pub const GSL_ROOT4_FLT_MIN = @as(f64, 3.2927225399135965e-10);
pub const GSL_ROOT5_FLT_MIN = @as(f64, 2.5944428542140822e-08);
pub const GSL_ROOT6_FLT_MIN = @as(f64, 4.7683715820312542e-07);
pub const GSL_LOG_FLT_MIN = -@as(f64, 8.7336544750553102e+01);
pub const GSL_FLT_MAX = @as(f64, 3.4028234663852886e+38);
pub const GSL_SQRT_FLT_MAX = @as(f64, 1.8446743523953730e+19);
pub const GSL_ROOT3_FLT_MAX = @as(f64, 6.9814635196223242e+12);
pub const GSL_ROOT4_FLT_MAX = @as(f64, 4.2949672319999986e+09);
pub const GSL_ROOT5_FLT_MAX = @as(f64, 5.0859007855960041e+07);
pub const GSL_ROOT6_FLT_MAX = @as(f64, 2.6422459233807749e+06);
pub const GSL_LOG_FLT_MAX = @as(f64, 8.8722839052068352e+01);
pub const GSL_SFLT_EPSILON = @as(f64, 4.8828125000000000e-04);
pub const GSL_SQRT_SFLT_EPSILON = @as(f64, 2.2097086912079612e-02);
pub const GSL_ROOT3_SFLT_EPSILON = @as(f64, 7.8745065618429588e-02);
pub const GSL_ROOT4_SFLT_EPSILON = @as(f64, 1.4865088937534013e-01);
pub const GSL_ROOT5_SFLT_EPSILON = @as(f64, 2.1763764082403100e-01);
pub const GSL_ROOT6_SFLT_EPSILON = @as(f64, 2.8061551207734325e-01);
pub const GSL_LOG_SFLT_EPSILON = -@as(f64, 7.6246189861593985e+00);
pub const GSL_MACH_EPS = GSL_DBL_EPSILON;
pub const GSL_SQRT_MACH_EPS = @as(f64, 3.2e-08);
pub const GSL_ROOT3_MACH_EPS = @as(f64, 1.0e-05);
pub const GSL_ROOT4_MACH_EPS = @as(f64, 0.000178);
pub const GSL_ROOT5_MACH_EPS = @as(f64, 0.00100);
pub const GSL_ROOT6_MACH_EPS = @as(f64, 0.00316);
pub const GSL_LOG_MACH_EPS = -@as(f64, 34.54);
pub const __GSL_PRECISION_H__ = "";
pub const _GSL_PREC_T_NUM = @as(c_int, 3);
pub const __GSL_NAN_H__ = "";
pub const GSL_POSINF = INFINITY;
pub const GSL_NEGINF = -INFINITY;
pub const GSL_NAN = NAN;
pub const GSL_POSZERO = @as(f64, 0.0);
pub const GSL_NEGZERO = -@as(f64, 0.0);
pub const __GSL_POW_INT_H__ = "";
pub const __GSL_MINMAX_H__ = "";
pub inline fn GSL_MAX(a: anytype, b: anytype) @TypeOf(if (a > b) a else b) {
    _ = &a;
    _ = &b;
    return if (a > b) a else b;
}
pub inline fn GSL_MIN(a: anytype, b: anytype) @TypeOf(if (a < b) a else b) {
    _ = &a;
    _ = &b;
    return if (a < b) a else b;
}
pub inline fn GSL_MAX_INT(a: anytype, b: anytype) @TypeOf(GSL_MAX(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MAX(a, b);
}
pub inline fn GSL_MIN_INT(a: anytype, b: anytype) @TypeOf(GSL_MIN(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MIN(a, b);
}
pub inline fn GSL_MAX_DBL(a: anytype, b: anytype) @TypeOf(GSL_MAX(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MAX(a, b);
}
pub inline fn GSL_MIN_DBL(a: anytype, b: anytype) @TypeOf(GSL_MIN(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MIN(a, b);
}
pub inline fn GSL_MAX_LDBL(a: anytype, b: anytype) @TypeOf(GSL_MAX(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MAX(a, b);
}
pub inline fn GSL_MIN_LDBL(a: anytype, b: anytype) @TypeOf(GSL_MIN(a, b)) {
    _ = &a;
    _ = &b;
    return GSL_MIN(a, b);
}
pub const M_SQRT3 = @as(f64, 1.73205080756887729352744634151);
pub const M_SQRTPI = @as(f64, 1.77245385090551602729816748334);
pub const M_LNPI = @as(f64, 1.14472988584940017414342735135);
pub const M_EULER = @as(f64, 0.57721566490153286060651209008);
pub inline fn GSL_IS_ODD(n: anytype) @TypeOf(n & @as(c_int, 1)) {
    _ = &n;
    return n & @as(c_int, 1);
}
pub inline fn GSL_IS_EVEN(n: anytype) @TypeOf(!(GSL_IS_ODD(n) != 0)) {
    _ = &n;
    return !(GSL_IS_ODD(n) != 0);
}
pub inline fn GSL_SIGN(x: anytype) @TypeOf(if (x >= @as(f64, 0.0)) @as(c_int, 1) else -@as(c_int, 1)) {
    _ = &x;
    return if (x >= @as(f64, 0.0)) @as(c_int, 1) else -@as(c_int, 1);
}
pub inline fn GSL_IS_REAL(x: anytype) @TypeOf(gsl_finite(x)) {
    _ = &x;
    return gsl_finite(x);
}
pub inline fn GSL_FN_EVAL(F: anytype, x: anytype) @TypeOf(F.*.function.*(x, F.*.params)) {
    _ = &F;
    _ = &x;
    return F.*.function.*(x, F.*.params);
}
pub inline fn GSL_FN_FDF_EVAL_F(FDF: anytype, x: anytype) @TypeOf(FDF.*.f.*(x, FDF.*.params)) {
    _ = &FDF;
    _ = &x;
    return FDF.*.f.*(x, FDF.*.params);
}
pub inline fn GSL_FN_FDF_EVAL_DF(FDF: anytype, x: anytype) @TypeOf(FDF.*.df.*(x, FDF.*.params)) {
    _ = &FDF;
    _ = &x;
    return FDF.*.df.*(x, FDF.*.params);
}
pub inline fn GSL_FN_FDF_EVAL_F_DF(FDF: anytype, x: anytype, y: anytype, dy: anytype) @TypeOf(FDF.*.fdf.*(x, FDF.*.params, y, dy)) {
    _ = &FDF;
    _ = &x;
    _ = &y;
    _ = &dy;
    return FDF.*.fdf.*(x, FDF.*.params, y, dy);
}
pub inline fn GSL_FN_VEC_EVAL(F: anytype, x: anytype, y: anytype) @TypeOf(F.*.function.*(x, y, F.*.params)) {
    _ = &F;
    _ = &x;
    _ = &y;
    return F.*.function.*(x, y, F.*.params);
}
pub const __GSL_BLAS_H__ = "";
pub const __GSL_COMPLEX_MATH_H__ = "";
pub const GSL_COMPLEX_ONE = gsl_complex_rect(@as(f64, 1.0), @as(f64, 0.0));
pub const GSL_COMPLEX_ZERO = gsl_complex_rect(@as(f64, 0.0), @as(f64, 0.0));
pub const GSL_COMPLEX_NEGONE = gsl_complex_rect(-@as(f64, 1.0), @as(f64, 0.0));
pub const __GSL_INTEGRATION_H__ = "";
pub const timeval = struct_timeval;
pub const timespec = struct_timespec;
pub const __pthread_internal_list = struct___pthread_internal_list;
pub const __pthread_internal_slist = struct___pthread_internal_slist;
pub const __pthread_mutex_s = struct___pthread_mutex_s;
pub const __pthread_rwlock_arch_t = struct___pthread_rwlock_arch_t;
pub const __pthread_cond_s = struct___pthread_cond_s;
pub const random_data = struct_random_data;
pub const drand48_data = struct_drand48_data;
pub const _G_fpos_t = struct__G_fpos_t;
pub const _G_fpos64_t = struct__G_fpos64_t;
pub const _IO_marker = struct__IO_marker;
pub const _IO_FILE = struct__IO_FILE;
pub const _IO_codecvt = struct__IO_codecvt;
pub const _IO_wide_data = struct__IO_wide_data;
pub const _IO_cookie_io_functions_t = struct__IO_cookie_io_functions_t;
pub const gsl_permutation_struct = struct_gsl_permutation_struct;
pub const gsl_block_long_double_struct = struct_gsl_block_long_double_struct;
pub const gsl_block_complex_long_double_struct = struct_gsl_block_complex_long_double_struct;
pub const gsl_block_struct = struct_gsl_block_struct;
pub const gsl_block_complex_struct = struct_gsl_block_complex_struct;
pub const gsl_block_float_struct = struct_gsl_block_float_struct;
pub const gsl_block_complex_float_struct = struct_gsl_block_complex_float_struct;
pub const gsl_block_ulong_struct = struct_gsl_block_ulong_struct;
pub const gsl_block_long_struct = struct_gsl_block_long_struct;
pub const gsl_block_uint_struct = struct_gsl_block_uint_struct;
pub const gsl_block_int_struct = struct_gsl_block_int_struct;
pub const gsl_block_ushort_struct = struct_gsl_block_ushort_struct;
pub const gsl_block_short_struct = struct_gsl_block_short_struct;
pub const gsl_block_uchar_struct = struct_gsl_block_uchar_struct;
pub const gsl_block_char_struct = struct_gsl_block_char_struct;
pub const CBLAS_ORDER = enum_CBLAS_ORDER;
pub const CBLAS_TRANSPOSE = enum_CBLAS_TRANSPOSE;
pub const CBLAS_UPLO = enum_CBLAS_UPLO;
pub const CBLAS_DIAG = enum_CBLAS_DIAG;
pub const CBLAS_SIDE = enum_CBLAS_SIDE;
pub const gsl_function_struct = struct_gsl_function_struct;
pub const gsl_function_fdf_struct = struct_gsl_function_fdf_struct;
pub const gsl_function_vec_struct = struct_gsl_function_vec_struct;
pub const gsl_integration_qawo_enum = enum_gsl_integration_qawo_enum;
