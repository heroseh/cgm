#ifndef CGM_H
#include "cgm.h"
#endif

// ===========================================================================
//
//
// Float Helpers
//
//
// ===========================================================================

float cgm_min(float a, float b) { return a < b ? a : b; }
float cgm_max(float a, float b) { return a > b ? a : b; }
float cgm_clamp(float v, float min, float max) { return v > max ? max : (v >= min ? v : min); }
float cgm_lerp(float from, float to, float t) { return (to - from) * t + from; }
float cgm_lerp_inv(float from, float to, float value) { return (value - from) / (to - from); }
float cgm_remap(float from_value, float from_min, float from_max, float to_min, float to_max) {
	float t = cgm_lerp_inv(from_min, from_max, from_value);
	return cgm_lerp(to_min, to_max, t);
}
CgmBool cgm_approx_eq(float a, float b) { return fabs(a - b) <= cgm_epsilon; }

float cgm_sign(float v) {
	if (v >= 0.f) return 1.f;
	return -1.f;
}

float cgm_round_to_multiple(float v, float multiple) {
	v += multiple * 0.5f;
	float rem = fmodf(v, multiple);
	if (v > 0.0f) {
		return v - rem;
	} else {
		return v - rem - multiple;
	}
}

float cgm_round_up_to_multiple(float v, float multiple) {
	float rem = fmodf(v, multiple);
	if (rem == 0.0) return v;
	if (v > 0.0) {
		return v + multiple - rem;
	} else {
		return v - rem;
	}
}

float cgm_round_down_to_multiple(float v, float multiple) {
	float rem = fmodf(v, multiple);
	if (rem == 0.0) return v;
	if (v > 0.0) {
		return v - rem;
	} else {
		return v - rem - multiple;
	}
}

typedef union _cgm_float_uint _cgm_float_uint;
union _cgm_float_uint {
	float f;
	uint32_t u;
};

float CgmF16_to_float(CgmF16 v) {
	if ((v.bits & 0x7c00) == 0x7c00) { // inf, -inf or nan
		if (v.bits & 0x03ff) return NAN;
		else if (v.bits & 0x8000) return -INFINITY;
		else return INFINITY;
	}

	_cgm_float_uint t1;
	uint32_t t2;
	uint32_t t3;

	t1.u = v.bits & 0x7fff;      // non-sign bits
	t2 = v.bits & 0x8000;        // sign bit
	t3 = v.bits & 0x7c00;        // exponent

	t1.u <<= 13;                 // align mantissa on MSB
	t2 <<= 16;                   // shift sign bit into position

	t1.u += 0x38000000;          // adjust bias

	t1.u = (t3 == 0 ? 0 : t1.u); // denormals-as-zero

	t1.u |= t2;                  // re-insert sign bit

	return t1.f;
}

CgmF16 CgmF16_from_float(float v) {
	if (isinf(v)) return (CgmF16) { .bits = v < 0.0 ? 0xfc00 : 0x7c00 };
	if (isnan(v)) return (CgmF16) { .bits = 0xffff };

	_cgm_float_uint vu = { .f = v };
	uint32_t t1;
	uint32_t t2;
	uint32_t t3;

	t1 = vu.u & 0x7fffffff;                // non-sign bits
	t2 = vu.u & 0x80000000;                // sign bit
	t3 = vu.u & 0x7f800000;                // exponent

	t1 >>= 13;                             // align mantissa on MSB
	t2 >>= 16;                             // shift sign bit into position

	t1 -= 0x1c000;                         // adjust bias

	t1 = (t3 < 0x38800000) ? 0 : t1;       // flush-to-zero
	t1 = (t3 > 0x8e000000) ? 0x7bff : t1;  // clamp-to-max
	t1 = (t3 == 0 ? 0 : t1);               // denormals-as-zero

	t1 |= t2;                              // re-insert sign bit

	return (CgmF16) { .bits = t1 };
}

CgmBool CgmF16_is_nan(CgmF16 v) {
	return (v.bits & 0x7c00) == 0x7c00 && v.bits & 0x03ff;
}

CgmBool CgmF16_is_inf(CgmF16 v) {
	return (v.bits & 0x7c00) == 0x7c00 && (v.bits & 0x03ff) == 0;

}

// ===========================================================================
//
//
// Vectors
//
//
// ===========================================================================

CgmBool CgmVec2_eq(CgmVec2 a, CgmVec2 b) { return a.x == b.x && a.y == b.y; }
CgmVec2 CgmVec2_copysign(CgmVec2 v, CgmVec2 copy) { return CgmVec2_init(copysignf(v.x, copy.x), copysignf(v.y, copy.y)); }

CgmVec2 CgmVec2_add(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(a.x + b.x, a.y + b.y); }
CgmVec2 CgmVec2_sub(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(a.x - b.x, a.y - b.y); }
CgmVec2 CgmVec2_mul(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(a.x * b.x, a.y * b.y); }
CgmVec2 CgmVec2_div(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(a.x / b.x, a.y / b.y); }
CgmVec2 CgmVec2_add_scalar(CgmVec2 v, float by) { return CgmVec2_init(v.x + by, v.y + by); }
CgmVec2 CgmVec2_sub_scalar(CgmVec2 v, float by) { return CgmVec2_init(v.x - by, v.y - by); }
CgmVec2 CgmVec2_mul_scalar(CgmVec2 v, float by) { return CgmVec2_init(v.x * by, v.y * by); }
CgmVec2 CgmVec2_div_scalar(CgmVec2 v, float by) { return CgmVec2_init(v.x / by, v.y / by); }

CgmVec2 CgmVec2_neg(CgmVec2 v) { return CgmVec2_init(-v.x, -v.y); }
float CgmVec2_len(CgmVec2 v) { return sqrtf((v.x * v.x) + (v.y * v.y)); }
CgmVec2 CgmVec2_norm(CgmVec2 v) {
	if (v.x == 0 && v.y == 0) return v;
	float k = 1.0 / sqrtf((v.x * v.x) + (v.y * v.y));
	return CgmVec2_init(v.x * k, v.y * k);
}

float CgmVec2_dot(CgmVec2 a, CgmVec2 b) {
	float p = 0.0;
	p += a.x * b.x;
	p += a.y * b.y;
	return p;
}

float CgmVec2_angle(CgmVec2 v) {
	return atan2f(-v.y, v.x);
}

CgmVec2 CgmVec2_mul_cross_scalar(CgmVec2 v, float s) {
	return CgmVec2_init(v.y * s, v.x * s);
}

float CgmVec2_mul_cross_vec(CgmVec2 a, CgmVec2 b) {
	return (a.x * b.y) - (a.y * b.x);
}

CgmVec2 CgmVec2_reflect(CgmVec2 v, CgmVec2 normal) {
	return CgmVec2_sub(v, CgmVec2_mul_scalar(normal, 2.0 * CgmVec2_dot(v, normal)));
}

CgmVec2 CgmVec2_absorb(CgmVec2 v, CgmVec2 normal) {
	CgmVec2 d = CgmVec2_mul_scalar(normal, CgmVec2_dot(v, normal));
	return CgmVec2_sub(v, d);
}

CgmVec2 CgmVec2_rotate(CgmVec2 v, float angle) {
	float co = cosf(angle);
	float si = sinf(angle);
	return CgmVec2_init(
		(v.x * si) - (v.y * co),
		(v.x * co) + (v.y * si)
	);
}

CgmVec2 CgmVec2_perp_left(CgmVec2 v) { return CgmVec2_init(v.y, -v.x); }
CgmVec2 CgmVec2_perp_right(CgmVec2 v) { return CgmVec2_init(-v.y, v.x); }

CgmVec2 CgmVec2_min(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(cgm_min(a.x, b.x), cgm_min(a.y, b.y)); }
CgmVec2 CgmVec2_max(CgmVec2 a, CgmVec2 b) { return CgmVec2_init(cgm_max(a.x, b.x), cgm_max(a.y, b.y)); }
CgmVec2 CgmVec2_clamp(CgmVec2 v, CgmVec2 min, CgmVec2 max) {
	return CgmVec2_init(cgm_clamp(v.x, min.x, max.x), cgm_clamp(v.y, min.y, max.y));
}
CgmVec2 CgmVec2_lerp(CgmVec2 from, CgmVec2 to, CgmVec2 t) {
	return CgmVec2_init(cgm_lerp(from.x, to.x, t.x), cgm_lerp(from.y, to.y, t.y));
}
CgmVec2 CgmVec2_sign(CgmVec2 v) {
	return CgmVec2_init(cgm_sign(v.x), cgm_sign(v.y));
}
CgmVec2 CgmVec2_round_to_multiple(CgmVec2 v, CgmVec2 multiple) {
	return CgmVec2_init(cgm_round_to_multiple(v.x, multiple.x), cgm_round_to_multiple(v.y, multiple.y));
}
CgmVec2 CgmVec2_round_up_to_multiple(CgmVec2 v, CgmVec2 multiple) {
	return CgmVec2_init(cgm_round_up_to_multiple(v.x, multiple.x), cgm_round_up_to_multiple(v.y, multiple.y));
}
CgmVec2 CgmVec2_round_down_to_multiple(CgmVec2 v, CgmVec2 multiple) {
	return CgmVec2_init(cgm_round_down_to_multiple(v.x, multiple.x), cgm_round_down_to_multiple(v.y, multiple.y));
}
CgmVec2 CgmVec2_abs(CgmVec2 v) { return CgmVec2_init(fabs(v.x), fabs(v.y)); }
CgmVec2 CgmVec2_floor(CgmVec2 v) { return CgmVec2_init(floorf(v.x), floorf(v.y)); }
CgmVec2 CgmVec2_ceil(CgmVec2 v) { return CgmVec2_init(ceilf(v.x), ceilf(v.y)); }
CgmVec2 CgmVec2_round(CgmVec2 v) { return CgmVec2_init(roundf(v.x), roundf(v.y)); }
CgmBool CgmVec2_approx_eq(CgmVec2 a, CgmVec2 b) { return cgm_approx_eq(a.x, b.x) && cgm_approx_eq(a.y, b.y); }


CgmBool CgmVec3_eq(CgmVec3 a, CgmVec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
CgmVec3 CgmVec3_copysign(CgmVec3 v, CgmVec3 copy) { return CgmVec3_init(copysignf(v.x, copy.x), copysignf(v.y, copy.y), copysignf(v.z, copy.z)); }

CgmVec3 CgmVec3_add(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(a.x + b.x, a.y + b.y, a.z + b.z); }
CgmVec3 CgmVec3_sub(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(a.x - b.x, a.y - b.y, a.z - b.z); }
CgmVec3 CgmVec3_mul(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(a.x * b.x, a.y * b.y, a.z * b.z); }
CgmVec3 CgmVec3_div(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(a.x / b.x, a.y / b.y, a.z / b.z); }
CgmVec3 CgmVec3_add_scalar(CgmVec3 v, float by) { return CgmVec3_init(v.x + by, v.y + by, v.z + by); }
CgmVec3 CgmVec3_sub_scalar(CgmVec3 v, float by) { return CgmVec3_init(v.x - by, v.y - by, v.z - by); }
CgmVec3 CgmVec3_mul_scalar(CgmVec3 v, float by) { return CgmVec3_init(v.x * by, v.y * by, v.z * by); }
CgmVec3 CgmVec3_div_scalar(CgmVec3 v, float by) { return CgmVec3_init(v.x / by, v.y / by, v.z / by); }

CgmVec3 CgmVec3_neg(CgmVec3 v) { return CgmVec3_init(-v.x, -v.y, -v.z); }
float CgmVec3_len(CgmVec3 v) { return sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z)); }

CgmVec3 CgmVec3_norm(CgmVec3 v) {
	float k = 1.0 / sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
	return CgmVec3_init(v.x * k, v.y * k, v.z * k);
}
float CgmVec3_dot(CgmVec3 a, CgmVec3 b) {
	float p = 0.0;
	p += a.x * b.x;
	p += a.y * b.y;
	p += a.z * b.z;
	return p;
}
CgmVec3 CgmVec3_mul_cross(CgmVec3 a, CgmVec3 b) {
	return CgmVec3_init(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

CgmVec3 CgmVec3_reflect(CgmVec3 v, CgmVec3 normal) {
	return CgmVec3_sub(v, CgmVec3_mul_scalar(normal, 2.0 * CgmVec3_dot(v, normal)));
}

CgmVec3 CgmVec3_absorb(CgmVec3 v, CgmVec3 normal) {
	CgmVec3 d = CgmVec3_mul_scalar(normal, CgmVec3_dot(v, normal));
	return CgmVec3_sub(v, d);
}

CgmVec3 CgmVec3_rotate(CgmVec3 point, CgmVec3 rotation_vector, float angle) {
	CgmVec3 rotation_vector_norm = CgmVec3_norm(rotation_vector);
	float cos_angle = cosf(angle);
	float dot = CgmVec3_dot(rotation_vector_norm, point);
	CgmVec3 cross = CgmVec3_mul_cross(rotation_vector_norm, point);
	float comp = (1.0f - cos_angle) * dot;
	CgmVec3 vec =
		CgmVec3_add(
			CgmVec3_add(
				CgmVec3_mul_scalar(rotation_vector_norm, comp),
				CgmVec3_mul_scalar(point, cos_angle)
			),
			CgmVec3_mul_scalar(cross, sinf(angle))
		);
	return vec;
}

CgmVec3 CgmVec3_clamp_magnitude(CgmVec3 v, float min, float max) {
	float len = CgmVec3_len(v);
	if (len == 0.f) return CgmVec3_zero;
	CgmVec3 norm = CgmVec3_mul_scalar(v, 1.f / len);
	return CgmVec3_mul_scalar(norm, cgm_clamp(len, min, max));
}


CgmVec3 CgmVec3_perp_left(CgmVec3 v) { return CgmVec3_init(-v.z, v.y, v.x); }
CgmVec3 CgmVec3_perp_right(CgmVec3 v) { return CgmVec3_init(v.z, v.y, -v.x); }
CgmVec3 CgmVec3_perp_forward(CgmVec3 v) { return CgmVec3_init(v.x, v.z, -v.y); }
CgmVec3 CgmVec3_perp_backward(CgmVec3 v) { return CgmVec3_init(v.x, -v.z, v.y); }

CgmVec3 CgmVec3_min(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(cgm_min(a.x, b.x), cgm_min(a.y, b.y), cgm_min(a.z, b.z)); }
CgmVec3 CgmVec3_max(CgmVec3 a, CgmVec3 b) { return CgmVec3_init(cgm_max(a.x, b.x), cgm_max(a.y, b.y), cgm_max(a.z, b.z)); }
CgmVec3 CgmVec3_clamp(CgmVec3 v, CgmVec3 min, CgmVec3 max) {
	return CgmVec3_init(cgm_clamp(v.x, min.x, max.x), cgm_clamp(v.y, min.y, max.y), cgm_clamp(v.z, min.z, max.z));
}
CgmVec3 CgmVec3_lerp(CgmVec3 from, CgmVec3 to, CgmVec3 t) {
	return CgmVec3_init(cgm_lerp(from.x, to.x, t.x), cgm_lerp(from.y, to.y, t.y), cgm_lerp(from.z, to.z, t.z));
}
CgmVec3 CgmVec3_sign(CgmVec3 v) {
	return CgmVec3_init(cgm_sign(v.x), cgm_sign(v.y), cgm_sign(v.z));
}
CgmVec3 CgmVec3_round_to_multiple(CgmVec3 v, CgmVec3 multiple) {
	return CgmVec3_init(cgm_round_to_multiple(v.x, multiple.x), cgm_round_to_multiple(v.y, multiple.y), cgm_round_to_multiple(v.z, multiple.z));
}
CgmVec3 CgmVec3_round_up_to_multiple(CgmVec3 v, CgmVec3 multiple) {
	return CgmVec3_init(cgm_round_up_to_multiple(v.x, multiple.x), cgm_round_up_to_multiple(v.y, multiple.y), cgm_round_up_to_multiple(v.z, multiple.z));
}
CgmVec3 CgmVec3_round_down_to_multiple(CgmVec3 v, CgmVec3 multiple) {
	return CgmVec3_init(cgm_round_down_to_multiple(v.x, multiple.x), cgm_round_down_to_multiple(v.y, multiple.y), cgm_round_down_to_multiple(v.z, multiple.z));
}
CgmVec3 CgmVec3_abs(CgmVec3 v) { return CgmVec3_init(fabs(v.x), fabs(v.y), fabs(v.z)); }
CgmVec3 CgmVec3_floor(CgmVec3 v) { return CgmVec3_init(floorf(v.x), floorf(v.y), floorf(v.z)); }
CgmVec3 CgmVec3_ceil(CgmVec3 v) { return CgmVec3_init(ceilf(v.x), ceilf(v.y), ceilf(v.z)); }
CgmVec3 CgmVec3_round(CgmVec3 v) { return CgmVec3_init(roundf(v.x), roundf(v.y), roundf(v.z)); }
CgmBool CgmVec3_approx_eq(CgmVec3 a, CgmVec3 b) { return cgm_approx_eq(a.x, b.x) && cgm_approx_eq(a.y, b.y) && cgm_approx_eq(a.z, b.z); }


CgmBool CgmVec4_eq(CgmVec4 a, CgmVec4 b) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }
CgmVec4 CgmVec4_copysign(CgmVec4 v, CgmVec4 copy) { return CgmVec4_init(copysignf(v.x, copy.x), copysignf(v.y, copy.y), copysignf(v.z, copy.z), copysignf(v.w, copy.w)); }

// ===========================================================================
//
//
// Matrices - row major order
//
//
// ===========================================================================

void CgmMat3x2_identity(CgmMat3x2* out) {
	out->row[0] = CgmVec2_init(1.0, 0.0);
	out->row[1] = CgmVec2_init(0.0, 1.0);
	out->row[2] = CgmVec2_init(0.0, 0.0);
}

void CgmMat3x2_identity_translate(CgmMat3x2* out, CgmVec2 v) {
	out->row[0] = CgmVec2_init(1.0, 0.0);
	out->row[1] = CgmVec2_init(0.0, 1.0);
	out->row[2] = v;
}

void CgmMat3x2_identity_scale(CgmMat3x2* out, CgmVec2 v) {
	out->row[0] = CgmVec2_init(v.x, 0.0);
	out->row[1] = CgmVec2_init(0.0, v.y);
	out->row[2] = CgmVec2_init(0.0, 0.0);
}

void CgmMat3x2_identity_rotate(CgmMat3x2* out, float angle) {
	float c = cosf(angle);
	float s = sinf(angle);
	out->row[0] = CgmVec2_init(c, -s);
	out->row[1] = CgmVec2_init(s, c);
	out->row[2] = CgmVec2_init(0.0, 0.0);
}

CgmVec2 CgmMat3x2_row(CgmMat3x2* m, uint32_t row_idx) {
	return m->row[row_idx];
}

CgmVec3 CgmMat3x2_column(CgmMat3x2* m, uint32_t column_idx) {
	return CgmVec3_init(
		m->row[0].a[column_idx],
		m->row[1].a[column_idx],
		m->row[2].a[column_idx]
	);
}

void CgmMat3x2_mul(CgmMat3x2* out, CgmMat3x2* a, CgmMat3x2* b) {
	out->row[0] = CgmVec2_init(
		(a->row[0].x * b->row[0].x) + (a->row[0].y * b->row[1].x),
		(a->row[0].x * b->row[0].y) + (a->row[0].y * b->row[1].y)
	);

	out->row[1] = CgmVec2_init(
		(a->row[1].x * b->row[0].x) + (a->row[1].y * b->row[1].x),
		(a->row[1].x * b->row[0].y) + (a->row[1].y * b->row[1].y)
	);

	out->row[2] = CgmVec2_init(
		(a->row[2].x * b->row[0].x) + (a->row[2].y * b->row[1].x) + b->row[2].x,
		(a->row[2].x * b->row[0].y) + (a->row[2].y * b->row[1].y) + b->row[2].y
	);
}

CgmVec2 CgmMat3x2_mul_point(CgmMat3x2* m, CgmVec2 pt) {
	return CgmVec2_init(
		(pt.x * m->row[0].x) + (pt.y * m->row[1].x) + m->row[2].x,
		(pt.x * m->row[0].y) + (pt.y * m->row[1].y) + m->row[2].y
	);
}

CgmVec2 CgmMat3x2_mul_vector(CgmMat3x2* m, CgmVec2 v) {
	return CgmVec2_init(
		(v.x * m->row[0].x) + (v.y * m->row[1].x),
		(v.x * m->row[0].y) + (v.y * m->row[1].y)
	);
}


void CgmMat4x4_identity(CgmMat4x4* out) {
	out->row[0] = CgmVec4_init(1.0, 0.0, 0.0, 0.0);
	out->row[1] = CgmVec4_init(0.0, 1.0, 0.0, 0.0);
	out->row[2] = CgmVec4_init(0.0, 0.0, 1.0, 0.0);
	out->row[3] = CgmVec4_init(0.0, 0.0, 0.0, 1.0);
}

void CgmMat4x4_identity_scale(CgmMat4x4* out, CgmVec3 v) {
	CgmMat4x4_identity(out);
	out->row[0].x = v.x;
	out->row[1].y = v.y;
	out->row[2].z = v.z;
}

void CgmMat4x4_identity_rotate(CgmMat4x4* out, CgmVec3 v, float angle) {
	float xx = v.x * v.x;
	float yy = v.y * v.y;
	float zz = v.z * v.z;

	float t = angle / 2.0;
	float ts = sinf(t);
	float tc = cosf(t);
	float sc = ts * tc;
	float sq = ts * ts;

	out->row[0] = CgmVec4_init(
		1.0 - (2.0 * (yy + zz) * sq),
		2.0 * ((v.x * v.y * sq) - (v.z * sc)),
		2.0 * ((v.x * v.z * sq) + (v.y * sc)),
		0.0
	);

	out->row[1] = CgmVec4_init(
		2.0 * ((v.x * v.y * sq) + (v.z * sc)),
		1.0 - (2.0 * (xx + zz) * sq),
		2.0 * ((v.y * v.z * sq) - (v.x * sc)),
		0.0
	);

	out->row[2] = CgmVec4_init(
		2.0 * ((v.x * v.z * sq) - (v.y * sc)),
		2.0 * ((v.y * v.z * sq) + (v.x * sc)),
		1.0 - (2.0 * (xx + yy) * sq),
		0.0
	);

	out->row[3] = CgmVec4_init(0.0, 0.0, 0.0, 1.0);
}

void CgmMat4x4_identity_translate(CgmMat4x4* out, CgmVec3 v) {
	CgmMat4x4_identity(out);
	out->row[3] = CgmVec4_init(v.x, v.y, v.z, 1.0);
}

void CgmMat4x4_scale(CgmMat4x4* m, CgmVec3 v) {
	CgmMat4x4 mat;
	CgmMat4x4_identity_scale(&mat, v);
	CgmMat4x4_mul(m, m, &mat);
}

void CgmMat4x4_rotate(CgmMat4x4* m, CgmVec3 v, float angle) {
	CgmMat4x4 mat;
	CgmMat4x4_identity_rotate(&mat, v, angle);
	CgmMat4x4_mul(m, m, &mat);
}

void CgmMat4x4_translate(CgmMat4x4* m, CgmVec3 v) {
	CgmMat4x4 mat;
	CgmMat4x4_identity_translate(&mat, v);
	CgmMat4x4_mul(m, m, &mat);
}

void CgmMat4x4_from_3x2(CgmMat4x4* out, CgmMat3x2* m) {
	out->row[0] = CgmVec4_init(m->row[0].x, m->row[0].y, 0.0, 0.0);
	out->row[1] = CgmVec4_init(m->row[1].x, m->row[1].y, 0.0, 0.0);
	out->row[2] = CgmVec4_init(0.0, 0.0, 1.0, 0.0);
	out->row[3] = CgmVec4_init(m->row[3].x, m->row[3].y, 0.0, 1.0);
}

CgmVec4 CgmMat4x4_row(CgmMat4x4* m, uint32_t row_idx) {
	return m->row[row_idx];
}

CgmVec4 CgmMat4x4_column(CgmMat4x4* m, uint32_t column_idx) {
	return CgmVec4_init(
		m->row[0].a[column_idx],
		m->row[1].a[column_idx],
		m->row[2].a[column_idx],
		m->row[3].a[column_idx]
	);
}

void CgmMat4x4_ortho(CgmMat4x4* out, float left, float right, float bottom, float top, float near, float far) {
	float diff_x = right - left;
	float diff_y = top - bottom;
	float diff_z = far - near;

	float tx = -((right + left) / diff_x);
	float ty = -((top + bottom) / diff_y);
	float tz = -((far + near) / diff_z);

	out->row[0] = CgmVec4_init(2.0 / diff_x, 0.0, 0.0, 0.0);
	out->row[1] = CgmVec4_init(0.0, 2.0 / diff_y, 0.0, 0.0);
	out->row[2] = CgmVec4_init(0.0, 0.0, -2.0 / diff_z, 0.0);
	out->row[3] = CgmVec4_init(tx, ty, tz, 1.0);
}

void CgmMat4x4_perspective(CgmMat4x4* out, float fovy, float aspect_ratio, float z_near, float z_far) {
	cgm_assert(aspect_ratio != 0.0, "aspect_ratio cannot be 0.0");
	cgm_assert(z_far != z_near, "z_near and z_far cannot be equal");

	float tan_half_fovy = tanf(fovy / 2.0);
	float a = 1.0 / tan_half_fovy;

	*out = (CgmMat4x4){0};
	out->row[0].x = a / aspect_ratio;
	out->row[1].y = a;
	out->row[2].z = -((z_far + z_near) / (z_far - z_near));
	out->row[2].w = -1.0;
	out->row[3].z = -((2.0 * z_far * z_near) / (z_far - z_near));
}

void CgmMat4x4_mul(CgmMat4x4* out, CgmMat4x4* a, CgmMat4x4* b) {
	out->row[0] = CgmVec4_init(
		a->row[0].x * b->row[0].x  +  a->row[0].y * b->row[1].x  +  a->row[0].z * b->row[2].x  +  a->row[0].w * b->row[3].x,
		a->row[0].x * b->row[0].y  +  a->row[0].y * b->row[1].y  +  a->row[0].z * b->row[2].y  +  a->row[0].w * b->row[3].y,
		a->row[0].x * b->row[0].z  +  a->row[0].y * b->row[1].z  +  a->row[0].z * b->row[2].z  +  a->row[0].w * b->row[3].z,
		a->row[0].x * b->row[0].w  +  a->row[0].y * b->row[1].w  +  a->row[0].z * b->row[2].w  +  a->row[0].w * b->row[3].w
	);

	out->row[1] = CgmVec4_init(
		a->row[1].x * b->row[0].x  +  a->row[1].y * b->row[1].x  +  a->row[1].z * b->row[2].x  +  a->row[1].w * b->row[3].x,
		a->row[1].x * b->row[0].y  +  a->row[1].y * b->row[1].y  +  a->row[1].z * b->row[2].y  +  a->row[1].w * b->row[3].y,
		a->row[1].x * b->row[0].z  +  a->row[1].y * b->row[1].z  +  a->row[1].z * b->row[2].z  +  a->row[1].w * b->row[3].z,
		a->row[1].x * b->row[0].w  +  a->row[1].y * b->row[1].w  +  a->row[1].z * b->row[2].w  +  a->row[1].w * b->row[3].w
	);

	out->row[2] = CgmVec4_init(
		a->row[2].x * b->row[0].x  +  a->row[2].y * b->row[1].x  +  a->row[2].z * b->row[2].x  +  a->row[2].w * b->row[3].x,
		a->row[2].x * b->row[0].y  +  a->row[2].y * b->row[1].y  +  a->row[2].z * b->row[2].y  +  a->row[2].w * b->row[3].y,
		a->row[2].x * b->row[0].z  +  a->row[2].y * b->row[1].z  +  a->row[2].z * b->row[2].z  +  a->row[2].w * b->row[3].z,
		a->row[2].x * b->row[0].w  +  a->row[2].y * b->row[1].w  +  a->row[2].z * b->row[2].w  +  a->row[2].w * b->row[3].w
	);

	out->row[3] = CgmVec4_init(
		a->row[3].x * b->row[0].x  +  a->row[3].y * b->row[1].x  +  a->row[3].z * b->row[2].x  +  a->row[3].w * b->row[3].x,
		a->row[3].x * b->row[0].y  +  a->row[3].y * b->row[1].y  +  a->row[3].z * b->row[2].y  +  a->row[3].w * b->row[3].y,
		a->row[3].x * b->row[0].z  +  a->row[3].y * b->row[1].z  +  a->row[3].z * b->row[2].z  +  a->row[3].w * b->row[3].z,
		a->row[3].x * b->row[0].w  +  a->row[3].y * b->row[1].w  +  a->row[3].z * b->row[2].w  +  a->row[3].w * b->row[3].w
	);
}

// ===========================================================================
//
//
// Collision Shapes 2D
//
//
// ===========================================================================

CgmVec2 CgmAabb2d_top_right(CgmAabb2d* a) {
	return CgmVec2_init(a->ex, a->y);
}

CgmVec2 CgmAabb2d_bottom_left(CgmAabb2d* a) {
	return CgmVec2_init(a->x, a->ey);
}

float CgmAabb2d_width(CgmAabb2d* a) {
	return a->ex - a->x;
}

float CgmAabb2d_height(CgmAabb2d* a) {
	return a->ey - a->y;
}

CgmVec2 CgmAabb2d_size(CgmAabb2d* a) {
	return CgmVec2_sub(a->max, a->min);
}

CgmVec2 CgmAabb2d_half_size(CgmAabb2d* a) {
	return CgmVec2_mul_scalar(CgmVec2_sub(a->max, a->min), 0.5);
}

CgmVec2 CgmAabb2d_center(CgmAabb2d* a) {
	return CgmVec2_add(a->min, CgmAabb2d_half_size(a));
}

// ===========================================================================
//
//
// Collision Checking 2D
//
//
// ===========================================================================

CgmBool cgm_2d_ray_vs_aabb(CgmRay2d* r, float r_max_distance, CgmAabb2d* a, CgmVec2* contact_normal_out, float* contact_distance_out) {
	float tx1 = (a->x - r->pos.x) / r->dir.x;
	float tx2 = (a->ex - r->pos.x) / r->dir.x;

	float tmin_x = cgm_min(tx1, tx2);
	float tmax = cgm_max(tx1, tx2);

	float ty1 = (a->y - r->pos.y) / r->dir.y;
	float ty2 = (a->ey - r->pos.y) / r->dir.y;

	float tmin_y = cgm_min(ty1, ty2);
	float tmin = cgm_max(tmin_x, tmin_y);
	tmax = cgm_min(tmax, cgm_max(ty1, ty2));

	if (tmax > 0.f && tmax >= tmin && tmin <= r_max_distance) {
		if (tmin_x > tmin_y) {
			*contact_normal_out = CgmVec2_init(-copysignf(1.0, r->dir.x), 0.0);
		} else {
			*contact_normal_out = CgmVec2_init(0.0, -copysignf(1.0, r->dir.y));
		}

		*contact_distance_out = tmin;
		return cgm_true;
	}

	return cgm_false;
}

CgmBool cgm_2d_swept_aabb_vs_aabb(CgmAabb2d* a, CgmVec2 a_dir, float a_max_distance, CgmAabb2d* b, CgmVec2* contact_normal_out, float* contact_distance_out) {
	if (a_dir.x == 0 && a_dir.y == 0) { return cgm_false; }

	CgmAabb2d tmp = *b;
	CgmVec2 a_half_size = CgmAabb2d_half_size(a);
	tmp.x -= a_half_size.x;
	tmp.y -= a_half_size.y;
	tmp.ex += a_half_size.x;
	tmp.ey += a_half_size.y;

	CgmRay2d ray = { .pos = CgmAabb2d_center(a), .dir = a_dir };
	return cgm_2d_ray_vs_aabb(&ray, a_max_distance, &tmp, contact_normal_out, contact_distance_out);
}

CgmBool cgm_2d_aabb_vs_swept_aabb(CgmAabb2d* a, CgmAabb2d* b, CgmVec2 b_dir, float b_max_distance, CgmVec2* contact_normal_out, float* contact_distance_out) {
	if (b_dir.x == 0 && b_dir.y == 0) { return cgm_false; }

	CgmAabb2d tmp = *a;
	CgmVec2 b_half_size = CgmAabb2d_half_size(b);
	tmp.x -= b_half_size.x;
	tmp.y -= b_half_size.y;
	tmp.ex += b_half_size.x;
	tmp.ey += b_half_size.y;

	CgmRay2d ray = { .pos = CgmAabb2d_center(b), .dir = b_dir };
	return cgm_2d_ray_vs_aabb(&ray, b_max_distance, &tmp, contact_normal_out, contact_distance_out);
}

CgmBool cgm_2d_swept_aabb_vs_swept_aabb(CgmAabb2d* a, CgmVec2 a_dir, float a_distance, CgmAabb2d* b, CgmVec2 b_dir, float b_distance, CgmVec2* contact_normal_out, float* contact_distance_out) {
	return cgm_2d_swept_aabb_vs_aabb(a, CgmVec2_sub(b_dir, a_dir), b_distance - a_distance, b, contact_normal_out, contact_distance_out);
}

CgmBool cgm_2d_aabb_vs_aabb(CgmAabb2d* a, CgmAabb2d* b, CgmVec2* contact_normal_out, float* contact_distance_out) {
	CgmAabb2d tmp = *b;
	CgmVec2 a_half_size = CgmAabb2d_half_size(a);
	tmp.x -= a_half_size.x;
	tmp.y -= a_half_size.y;
	tmp.ex += a_half_size.x;
	tmp.ey += a_half_size.y;

	return cgm_2d_aabb_vs_pt(&tmp, CgmAabb2d_center(a), contact_normal_out, contact_distance_out);
}

CgmBool cgm_2d_aabb_vs_pt(CgmAabb2d* a, CgmVec2 pt, CgmVec2* contact_normal_out, float* contact_distance_out) {
	CgmVec2 vec = CgmVec2_sub(pt, CgmAabb2d_center(a));
	CgmVec2 relative_pt = CgmVec2_sub(CgmAabb2d_half_size(a), CgmVec2_abs(vec));

	if (relative_pt.x < 0 || relative_pt.y < 0) {
		return cgm_false;
	}

	if (relative_pt.x < relative_pt.y) {
		if (contact_distance_out) *contact_distance_out = relative_pt.x;
		if (contact_normal_out) *contact_normal_out = CgmVec2_init(copysignf(1.0, vec.x), 0.0);
	} else {
		if (contact_distance_out) *contact_distance_out = relative_pt.y;
		if (contact_normal_out) *contact_normal_out = CgmVec2_init(0.0, copysignf(1.0, vec.y));
	}

	return cgm_true;
}

// ===========================================================================
//
//
// Collision Shapes 3D
//
//
// ===========================================================================

float CgmAabb3d_width(CgmAabb3d* a) { return a->ex - a->x; }
float CgmAabb3d_height(CgmAabb3d* a) { return a->ey - a->y; }
float CgmAabb3d_depth(CgmAabb3d* a) { return a->ez - a->z; }
CgmVec3 CgmAabb3d_size(CgmAabb3d* a) { return CgmVec3_init(a->ex - a->x, a->ey - a->y, a->ez - a->z); }
CgmVec3 CgmAabb3d_half_size(CgmAabb3d* a) { return CgmVec3_mul_scalar(CgmAabb3d_size(a), 0.5f); }
CgmVec3 CgmAabb3d_center(CgmAabb3d* a) { return CgmVec3_add(a->min, CgmVec3_mul_scalar(CgmAabb3d_size(a), 0.5f)); }
CgmVec3 CgmAabb3d_clamp_pt(CgmAabb3d* a, CgmVec3 pt) { return CgmVec3_clamp(pt, a->min, a->max); }
CgmVec3 CgmAabb3d_closest_pt_on_ray(CgmAabb3d* a, CgmVec3 origin, CgmVec3 direction, float min, float max) {
	CgmVec3 vertices[] = {
		CgmVec3_init(a->x,   a->y,  a->z),
		CgmVec3_init(a->ex,  a->y,  a->z),
		CgmVec3_init(a->x,  a->ey,  a->z),
		CgmVec3_init(a->x,   a->y, a->ez),
		CgmVec3_init(a->ex, a->ey,  a->z),
		CgmVec3_init(a->x,  a->ey, a->ez),
		CgmVec3_init(a->ex,  a->y, a->ez),
		CgmVec3_init(a->ex, a->ey, a->ez),
	};

	float min_distance_squared = INFINITY;
	CgmVec3 closest_pt;
	for (uint32_t i = 0; i < 8; i += 1) {
		CgmVec3 vertex = vertices[i];
		CgmVec3 vec_to_vertex = CgmVec3_sub(vertex, origin);
		float offset_to_perp_pt = CgmVec3_dot(vec_to_vertex, direction);
		offset_to_perp_pt = cgm_clamp(offset_to_perp_pt, min, max);
		CgmVec3 perp_pt = CgmVec3_add(origin, CgmVec3_mul_scalar(direction, offset_to_perp_pt));
		CgmVec3 perp_vec = CgmVec3_sub(vertex, perp_pt);
		float distance_squared = CgmVec3_dot(perp_vec, perp_vec);
		if (distance_squared < min_distance_squared) {
			closest_pt = perp_pt;
			min_distance_squared = distance_squared;
		}
	}

	return closest_pt;
}

void CgmSphere_aabb(CgmSphere* sphere, CgmAabb3d* aabb_in_out) {
	CgmVec3 min = CgmVec3_sub_scalar(sphere->center_pos, sphere->radius);
	CgmVec3 max = CgmVec3_add_scalar(sphere->center_pos, sphere->radius);
	aabb_in_out->min = min;
	aabb_in_out->max = max;
}

void CgmCapsule3d_aabb(CgmCapsule3d* capsule, CgmAabb3d* aabb_in_out) {
	CgmVec3 offset = CgmVec3_mul_scalar(capsule->direction, capsule->half_length);
	CgmVec3 start = CgmVec3_add(capsule->center_pos, offset);
	CgmVec3 end = CgmVec3_add(capsule->center_pos, CgmVec3_neg(offset));

	CgmSphere sphere_start = { .radius = capsule->radius, .center_pos = start };
	CgmSphere sphere_end = { .radius = capsule->radius, .center_pos = end };
	CgmSphere_aabb(&sphere_start, aabb_in_out);
	CgmSphere_aabb(&sphere_end, aabb_in_out);
}

CgmVec3 CgmCapsule3d_get_projected_pt(CgmCapsule3d* capsule, CgmVec3 pt) {
	CgmVec3 vec_to_pt = CgmVec3_sub(pt, capsule->center_pos);
	float offset = CgmVec3_dot(capsule->direction, vec_to_pt);
	offset = cgm_clamp(offset, -capsule->half_length, capsule->half_length);
	return CgmVec3_add(capsule->center_pos, CgmVec3_mul_scalar(capsule->direction, offset));
}

void CgmTriangle3d_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb_in_out) {
	aabb_in_out->x = cgm_min(triangle->a.x, cgm_min(triangle->b.x, triangle->c.x));
	aabb_in_out->ex = cgm_max(triangle->a.x, cgm_max(triangle->b.x, triangle->c.x));
	aabb_in_out->y = cgm_min(triangle->a.y, cgm_min(triangle->b.y, triangle->c.y));
	aabb_in_out->ey = cgm_max(triangle->a.y, cgm_max(triangle->b.y, triangle->c.y));
	aabb_in_out->z = cgm_min(triangle->a.z, cgm_min(triangle->b.z, triangle->c.z));
	aabb_in_out->ez = cgm_max(triangle->a.z, cgm_max(triangle->b.z, triangle->c.z));
}

void CgmTriangle3d_offset(CgmTriangle3d* triangle, CgmVec3 by) {
	triangle->a = CgmVec3_add(triangle->a, by);
	triangle->b = CgmVec3_add(triangle->b, by);
	triangle->c = CgmVec3_add(triangle->c, by);
}

CgmVec3 cgm_line_get_projected_pt(CgmVec3 start, CgmVec3 end, CgmVec3 pt) {
	CgmVec3 vec = CgmVec3_sub(end, start);
	float t = CgmVec3_dot(CgmVec3_sub(pt, start), vec) / CgmVec3_dot(vec, vec);
	return CgmVec3_add(start, CgmVec3_mul_scalar(vec, cgm_clamp(t, 0.f, 1.f)));
}


// ===========================================================================
//
//
// Collision Checking 3D
//
//
// ===========================================================================

CgmBool cgm_3d_pt_vs_aabb(CgmVec3 pt, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmVec3 vec = CgmVec3_sub(pt, CgmAabb3d_center(aabb));
	CgmVec3 relative_pt = CgmVec3_sub(CgmAabb3d_half_size(aabb), CgmVec3_abs(vec));

	if (relative_pt.x < 0 || relative_pt.y < 0 || relative_pt.z < 0) {
		return cgm_false;
	}

	if (relative_pt.x < relative_pt.y && relative_pt.x < relative_pt.z) {
		if (contact_distance_out) *contact_distance_out = relative_pt.x;
		if (contact_normal_out) *contact_normal_out = CgmVec3_init(-copysignf(1.0, vec.x), 0.f, 0.f);
	} else if (relative_pt.y < relative_pt.x && relative_pt.y < relative_pt.z) {
		if (contact_distance_out) *contact_distance_out = relative_pt.y;
		if (contact_normal_out) *contact_normal_out = CgmVec3_init(0.0, -copysignf(1.0, vec.y), 0.f);
	} else {
		if (contact_distance_out) *contact_distance_out = relative_pt.z;
		if (contact_normal_out) *contact_normal_out = CgmVec3_init(0.0, 0.f, -copysignf(1.0, vec.z));
	}

	return cgm_true;
}

CgmBool cgm_3d_pt_vs_sphere(CgmVec3 pt, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmVec3 vec = CgmVec3_sub(pt, sphere->center_pos);
	float vec_len_squared = CgmVec3_dot(vec, vec);
	float radius_squared = sphere->radius * sphere->radius;
	if (vec_len_squared <= radius_squared) {
		if (!contact_normal_out && !contact_distance_out)
			return cgm_true;

		float distance = sqrtf(vec_len_squared);
		if (contact_distance_out) *contact_distance_out = sphere->radius - distance;
		if (contact_normal_out) {
			if (distance == 0.f) *contact_normal_out = CgmVec3_init(0.f, 0.f, 1.f);
			else {
				// CgmVec3_norm(vec)
				float k = 1.0 / sqrtf(vec_len_squared);
				*contact_normal_out = CgmVec3_init(vec.x * k, vec.y * k, vec.z * k);
			}
		}
		return cgm_true;
	}

	return cgm_false;
}

CgmBool cgm_3d_pt_vs_capsule(CgmVec3 pt, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmVec3 vec_to_pt = CgmVec3_sub(pt, capsule->center_pos);
	float offset = CgmVec3_dot(capsule->direction, vec_to_pt);
	offset = cgm_clamp(offset, -capsule->half_length, capsule->half_length);

	CgmVec3 projected_pt = CgmVec3_add(capsule->center_pos, CgmVec3_mul_scalar(capsule->direction, offset));

	CgmSphere sphere = { .radius = capsule->radius, .center_pos = projected_pt };
	return cgm_3d_pt_vs_sphere(pt, &sphere, contact_normal_out, contact_distance_out);
}

CgmBool cgm_3d_aabb_vs_pt(CgmAabb3d* aabb, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_pt_vs_aabb(pt, aabb, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_aabb_vs_aabb(CgmAabb3d* a, CgmAabb3d* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmAabb3d tmp = *b;
	CgmVec3 a_half_size = CgmAabb3d_half_size(a);
	tmp.x -= a_half_size.x;
	tmp.y -= a_half_size.y;
	tmp.z -= a_half_size.z;
	tmp.ex += a_half_size.x;
	tmp.ey += a_half_size.y;
	tmp.ez += a_half_size.z;

	return cgm_3d_aabb_vs_pt(&tmp, CgmAabb3d_center(a), contact_normal_out, contact_distance_out);
}

CgmBool cgm_3d_aabb_vs_sphere(CgmAabb3d* aabb, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	if (cgm_3d_aabb_vs_pt(aabb, sphere->center_pos, NULL, NULL)) {
		//
		// the sphere center is inside the AABB if we do not want
		// any of the contact info returned we can return success early.
		//
		if (!contact_distance_out && !contact_normal_out)
			return cgm_true;

		//
		// determine the closest edge to get the normal and distance from that edge.
		//
		CgmVec3 vec_to_sphere = CgmVec3_sub(sphere->center_pos, CgmAabb3d_center(aabb));
		CgmVec3 vec_to_sphere_abs = CgmVec3_abs(vec_to_sphere);
		CgmVec3 half_size = CgmAabb3d_half_size(aabb);
		CgmVec3 diff = CgmVec3_sub(half_size, vec_to_sphere_abs);

		CgmVec3 normal;
		float distance_to_edge;
		if (diff.x < diff.y && diff.x < diff.z) {
			normal = CgmVec3_init(copysignf(1.f, -vec_to_sphere.x), 0.f, 0.f);
			distance_to_edge = diff.x;
		} else if (diff.y < diff.x && diff.y < diff.z) {
			normal = CgmVec3_init(0.f, copysignf(1.f, -vec_to_sphere.y), 0.f);
			distance_to_edge = diff.y;
		} else {
			normal = CgmVec3_init(0.f, 0.f, copysignf(1.f, -vec_to_sphere.z));
			distance_to_edge = diff.z;
		}

		if (contact_distance_out) *contact_distance_out = sphere->radius + distance_to_edge;
		if (contact_normal_out) *contact_normal_out = normal;
		return cgm_true;
	} else {
		CgmVec3 closest_pt = CgmAabb3d_clamp_pt(aabb, sphere->center_pos);
		CgmVec3 vec = CgmVec3_sub(sphere->center_pos, closest_pt);

		float distance_squared = CgmVec3_dot(vec, vec);
		float radius_squared = sphere->radius * sphere->radius;

		if (distance_squared <= radius_squared) {
			if (!contact_distance_out && !contact_normal_out)
				return cgm_true;

			float distance = sqrtf(distance_squared);
			if (contact_distance_out) *contact_distance_out = sphere->radius - distance;
			if (contact_normal_out) *contact_normal_out = CgmVec3_mul_scalar(vec, -(1.f / distance));
			return cgm_true;
		}
	}
	return cgm_false;
}

CgmBool cgm_3d_aabb_vs_capsule(CgmAabb3d* aabb, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmVec3 aabb_center = CgmAabb3d_center(aabb);
	CgmVec3 aabb_half_size = CgmAabb3d_half_size(aabb);

	// translate the capsule into relative space from the AABB center.
	// and the start and end in relative space as well.
	CgmVec3 capsule_rel_pos = CgmVec3_sub(capsule->center_pos, aabb_center);
	CgmVec3 capsule_offset = CgmVec3_mul_scalar(capsule->direction, capsule->half_length);
	CgmVec3 capsule_start = CgmVec3_add(capsule_rel_pos, capsule_offset);
	CgmVec3 capsule_end = CgmVec3_add(capsule_rel_pos, CgmVec3_neg(capsule_offset));

	// calculate the vector that goes down the capsule from the start to the end
	CgmVec3 capsule_edge_vector = CgmVec3_sub(capsule_end, capsule_start);

	// compute the face normals of the AABB, because the AABB
	// is at center, and of course axis aligned, we know that
	// it's normals are the X, Y and Z axis.
	CgmVec3 u0 = CgmVec3_init(1.0f, 0.0f, 0.0f);
	CgmVec3 u1 = CgmVec3_init(0.0f, 1.0f, 0.0f);
	CgmVec3 u2 = CgmVec3_init(0.0f, 0.0f, 1.0f);

	//
	// all combinations of the edges from both shapes paired together.
	// we will use the cross product on these pairs to get the
	// separation axis
	CgmVec3 edge_axises[][2] = {
		{ u0, capsule_edge_vector },
		{ u1, capsule_edge_vector },
		{ u2, capsule_edge_vector },
	};

	//
	// we need a separation axis to make sure the start and end of the capsule
	// with the nearest edges of the AABB are separated.
	CgmVec3 capsule_start_to_closest_aabb_pt_vec =
		CgmVec3_norm(CgmVec3_sub(capsule_start, CgmVec3_clamp(capsule_start, CgmVec3_neg(aabb_half_size), aabb_half_size)));
	CgmVec3 capsule_end_to_closest_aabb_pt_vec =
		CgmVec3_norm(CgmVec3_sub(capsule_end, CgmVec3_clamp(capsule_end, CgmVec3_neg(aabb_half_size), aabb_half_size)));

	//
	// all of the separation axis that do not need a cross product to work out.
	CgmVec3 face_axises[] = {
		u0,
		u1,
		u2,
		capsule_start_to_closest_aabb_pt_vec,
		capsule_end_to_closest_aabb_pt_vec,
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;

	for (uint32_t i = 0; i < 8; i += 1) {
		CgmVec3 axis;
		if (i < 3) {
			axis = edge_axises[i][0];
			CgmVec3 separation_axis = CgmVec3_mul_cross(axis, edge_axises[i][1]);
			float squared_distance = CgmVec3_dot(separation_axis, separation_axis);
			// if the two axis we just cross product together are parallel, we can skip these.
			// this is only because the edge axis is also use as a face axis -> u0, u1, u2
			// and will be tested later.
			if (cgm_approx_eq(squared_distance, 0.f)) continue;
			axis = CgmVec3_norm(separation_axis);
		} else {
			axis = face_axises[i - 3];
		}

		// project the capsule start and end onto the axis
		float p_start = CgmVec3_dot(capsule_start, axis);
		float p_end = CgmVec3_dot(capsule_end, axis);

		// project the half size of the AABB onto the separating axis
		// since we are working in relative space from the AABB center.
		float r = aabb_half_size.x * fabsf(CgmVec3_dot(u0, axis)) +
					aabb_half_size.y * fabsf(CgmVec3_dot(u1, axis)) +
					aabb_half_size.z * fabsf(CgmVec3_dot(u2, axis));

		// now do the actual test to see if either of
		// the most extreme of the capsule ends + their radius intersects r.
		float max = cgm_max(p_start, p_end) + capsule->radius;
		float min = cgm_min(p_start, p_end) - capsule->radius;
		float t = cgm_max(-max, min);
		if (t > r) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		float distance = r - t;
		if (distance < min_axis_distance) {
			//
			// we have found a axis with the new minimum separation needed to separate.
			// we need to flip the axis if the capsule is in the positive half of the AABB (AABB center is 0 AKA origin).
			min_axis_distance = distance;
			min_axis = min + ((max - min) * 0.5f) > 0.f ? CgmVec3_neg(axis) : axis;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_aabb_vs_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmVec3 aabb_center = CgmAabb3d_center(aabb);
	CgmVec3 aabb_half_size = CgmAabb3d_half_size(aabb);

	// translate triangle into relative space from the AABB center.
	CgmVec3 rp0 = CgmVec3_sub(triangle->a, aabb_center);
	CgmVec3 rp1 = CgmVec3_sub(triangle->b, aabb_center);
	CgmVec3 rp2 = CgmVec3_sub(triangle->c, aabb_center);

	// compute edge vectors for triangle
	CgmVec3 e0 = CgmVec3_sub(rp1, rp0);
	CgmVec3 e1 = CgmVec3_sub(rp2, rp1);
	CgmVec3 e2 = CgmVec3_sub(rp0, rp2);

	CgmVec3 triangle_face_normal = CgmVec3_norm(CgmVec3_mul_cross(e0, e1));

	// compute the face normals of the AABB, because the AABB
	// is at center, and of course axis aligned, we know that
	// it's normals are the X, Y and Z axis.
	CgmVec3 u0 = CgmVec3_init(1.0f, 0.0f, 0.0f);
	CgmVec3 u1 = CgmVec3_init(0.0f, 1.0f, 0.0f);
	CgmVec3 u2 = CgmVec3_init(0.0f, 0.0f, 1.0f);

	//
	// all combinations of the edges from both shapes paired together.
	// we will use the cross product on these pairs to get the
	// separation axis
	CgmVec3 edge_axises[][2] = {
		{ u0, e0 },
		{ u0, e1 },
		{ u0, e2 },
		{ u1, e0 },
		{ u1, e1 },
		{ u1, e2 },
		{ u2, e0 },
		{ u2, e1 },
		{ u2, e2 },
	};

	//
	// all of the separation axis that do not need a cross product to work out.
	CgmVec3 face_axises[] = {
		triangle_face_normal,
		u0,
		u1,
		u2,
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;

	for (uint32_t i = 0; i < 13; i += 1) {
		CgmVec3 axis;
		if (i < 9) {
			axis = edge_axises[i][0];
			CgmVec3 separation_axis = CgmVec3_mul_cross(axis, edge_axises[i][1]);
			float squared_distance = CgmVec3_dot(separation_axis, separation_axis);
			// if the two axis we just cross product together are parallel, we can skip these.
			// this is only because the edge axis is also use as a face axis -> u0, u1, u2
			// and will be tested later.
			if (cgm_approx_eq(squared_distance, 0.f)) continue;
			axis = CgmVec3_norm(separation_axis);
		} else {
			axis = face_axises[i - 9];
		}
		// project all 3 vertices of the triangle onto the separation axis
		float p0 = CgmVec3_dot(rp0, axis);
		float p1 = CgmVec3_dot(rp1, axis);
		float p2 = CgmVec3_dot(rp2, axis);

		// project the half size of the AABB onto the separating axis
		// since we are working in relative space from the AABB center.
		float r = aabb_half_size.x * fabsf(CgmVec3_dot(u0, axis)) +
					aabb_half_size.y * fabsf(CgmVec3_dot(u1, axis)) +
					aabb_half_size.z * fabsf(CgmVec3_dot(u2, axis));

		// now do the actual test to see if either of
		// the most extreme of the triangle points intersect with r.
		float max = cgm_max(cgm_max(p0, p1), p2);
		float min = cgm_min(cgm_min(p0, p1), p2);
		float t = cgm_max(-max, min);
		if (t > r) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		float distance = r - t;
		if (distance < min_axis_distance) {
			//
			// we have found a axis with the new minimum separation needed to separate.
			// we need to flip the axis if the triangle is in the positive half of the AABB (AABB center is 0 AKA origin).
			min_axis_distance = distance;
			min_axis = min + ((max - min) * 0.5f) > 0.f ? CgmVec3_neg(axis) : axis;;
		}
	}

	//
	// when using single side triangles, we only want to provide
	// a collision separation when the minimum separation axis
	// is the face normal of the triangle
	if (triangle_single_side) {
		float d = CgmVec3_dot(min_axis, triangle_face_normal);
		if (d < 0.70710678118) {
			return cgm_false;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_aabb_vs_single_side_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_aabb_vs_triangle(aabb, triangle, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_aabb_vs_double_side_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_aabb_vs_triangle(aabb, triangle, contact_normal_out, contact_distance_out, cgm_false);
}

static inline CgmBool _cgm_3d_shape_vs_mesh(void* shape, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool(*fn)()) {
	CgmTriangle3d triangle;
	CgmVec3 best_contact_normal;
	float max_contact_distance = -INFINITY;
	CgmVec3 contact_normal;
	float contact_distance;
	for (uint32_t i = 0; i < mesh->triangles_count; i += 1) {
		triangle = mesh->triangles[i];
		CgmTriangle3d_offset(&triangle, mesh->center_pos);
		if (fn(shape, &triangle, &contact_normal, &contact_distance, cgm_true)) {
			if (!contact_normal_out && !contact_distance_out) return cgm_true;

			if (contact_distance > max_contact_distance) {
				max_contact_distance = contact_distance;
				best_contact_normal = contact_normal;
			}
		}
	}

	if (isinf(max_contact_distance)) return cgm_false;

	if (contact_normal_out) *contact_normal_out = best_contact_normal;
	if (contact_distance_out) *contact_distance_out = contact_distance;
	return cgm_true;
}

CgmBool cgm_3d_aabb_vs_mesh(CgmAabb3d* aabb, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return _cgm_3d_shape_vs_mesh(aabb, mesh, contact_normal_out, contact_distance_out, (CgmBool(*)())cgm_3d_aabb_vs_triangle);
}

CgmBool cgm_3d_sphere_vs_pt(CgmSphere* sphere, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_pt_vs_sphere(pt, sphere, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_sphere_vs_sphere(CgmSphere* a, CgmSphere* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmSphere tmp = { .radius = a->radius + b->radius, .center_pos = b->center_pos };
	return cgm_3d_pt_vs_sphere(a->center_pos, &tmp, contact_normal_out, contact_distance_out);
}

CgmBool cgm_3d_sphere_vs_aabb(CgmSphere* sphere, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_aabb_vs_sphere(aabb, sphere, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_sphere_vs_capsule(CgmSphere* sphere, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmVec3 projected_pt = CgmCapsule3d_get_projected_pt(capsule, sphere->center_pos);
	CgmSphere tmp = { .radius = capsule->radius, .center_pos = projected_pt };
	return cgm_3d_sphere_vs_sphere(sphere, &tmp, contact_normal_out, contact_distance_out);
}

CgmBool cgm_3d_sphere_vs_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	// translate triangle into relative space from the AABB center.
	CgmVec3 rp0 = CgmVec3_sub(triangle->a, sphere->center_pos);
	CgmVec3 rp1 = CgmVec3_sub(triangle->b, sphere->center_pos);
	CgmVec3 rp2 = CgmVec3_sub(triangle->c, sphere->center_pos);

	// compute edge vectors for triangle
	CgmVec3 e0 = CgmVec3_sub(triangle->b, triangle->a);
	CgmVec3 e1 = CgmVec3_sub(triangle->c, triangle->b);
	CgmVec3 e2 = CgmVec3_sub(triangle->a, triangle->c);

	CgmVec3 triangle_face_normal = CgmVec3_norm(CgmVec3_mul_cross(e0, e1));

	//
	// all of the separation axis that do not need a cross product to work out.
	CgmVec3 face_axises[] = {
		triangle_face_normal,

		//
		// these three are essentially computing the closest
		// vector (the vector that is perpendicular) from each edge
		// to the sphere. CgmVec3_zero is used as the sphere center
		// position in relative space.
		CgmVec3_neg(cgm_line_get_projected_pt(rp0, rp1, CgmVec3_zero)),
		CgmVec3_neg(cgm_line_get_projected_pt(rp1, rp2, CgmVec3_zero)),
		CgmVec3_neg(cgm_line_get_projected_pt(rp2, rp0, CgmVec3_zero)),
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;

	for (uint32_t i = 0; i < 4; i += 1) {
		CgmVec3 axis = face_axises[i];
		float len_sq = CgmVec3_dot(axis, axis);
		if (len_sq == 0.f) continue;
		// normalize
		float k = 1.f / sqrtf(len_sq);
		axis = CgmVec3_mul_scalar(axis, k);

		// project all 3 vertices of the triangle onto the separation axis
		float p0 = CgmVec3_dot(rp0, axis);
		float p1 = CgmVec3_dot(rp1, axis);
		float p2 = CgmVec3_dot(rp2, axis);

		// now do the actual test to see if either of
		// the most extreme of the triangle points intersect with sphere->radius.
		float max = cgm_max(cgm_max(p0, p1), p2);
		float min = cgm_min(cgm_min(p0, p1), p2);
		float t = cgm_max(-max, min);
		if (t > sphere->radius) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		float distance = sphere->radius - t;
		if (distance < min_axis_distance) {
			//
			// we have found a axis with the new minimum separation needed to separate.
			// we need to flip the axis if the triangle is in the positive half of the sphere (sphere center is 0 AKA origin).
			min_axis_distance = distance;
			min_axis = min + ((max - min) * 0.5f) > 0.f ? CgmVec3_neg(axis) : axis;
		}
	}

	//
	// when using single side triangles, we only want to provide
	// a collision separation when the minimum separation axis
	// is the face normal of the triangle
	if (triangle_single_side) {
		float d = CgmVec3_dot(min_axis, triangle_face_normal);
		if (!cgm_approx_eq(d, 1.f)) {
			return cgm_false;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_sphere_vs_single_side_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_sphere_vs_triangle(sphere, triangle, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_sphere_vs_double_side_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_sphere_vs_triangle(sphere, triangle, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_sphere_vs_mesh(CgmSphere* sphere, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return _cgm_3d_shape_vs_mesh(sphere, mesh, contact_normal_out, contact_distance_out, (CgmBool(*)())cgm_3d_sphere_vs_triangle);
}

CgmBool cgm_3d_capsule_vs_pt(CgmCapsule3d* capsule, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_pt_vs_capsule(pt, capsule, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_capsule_vs_capsule(CgmCapsule3d* a, CgmCapsule3d* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	//
	// calculate the start and end positions of capsule A and capsule B.
	// as well as an edge vector that goes from the start to end of each capsule.
	CgmVec3 a_offset = CgmVec3_mul_scalar(a->direction, a->half_length);
	CgmVec3 a_start = CgmVec3_add(a->center_pos, a_offset);
	CgmVec3 a_end = CgmVec3_add(a->center_pos, CgmVec3_neg(a_offset));
	CgmVec3 a_edge_vector = CgmVec3_sub(a_end, a_start);

	CgmVec3 b_offset = CgmVec3_mul_scalar(b->direction, b->half_length);
	CgmVec3 b_start = CgmVec3_add(b->center_pos, b_offset);
	CgmVec3 b_end = CgmVec3_add(b->center_pos, CgmVec3_neg(b_offset));
	CgmVec3 b_edge_vector = CgmVec3_sub(b_end, b_start);

	//
	// all combinations of the edges from both shapes paired together.
	// we will use the cross product on these pairs to get the
	// separation axis
	CgmVec3 edge_axises[][2] = {
		{ a_edge_vector, b_edge_vector },
	};

	//
	// for each start and end of a capsule, project that point
	// on to the other capsule and make a vector towards that
	// point. then find the closest vector so we can use this
	// as a separation axis
	CgmVec3 closest_end_to_projected_pt_vec;
	{
		float closest_distance_squared = INFINITY;
		CgmVec3 end_to_projected_pt_vecs[] = {
			CgmVec3_sub(CgmCapsule3d_get_projected_pt(a, b_start), b_start),
			CgmVec3_sub(CgmCapsule3d_get_projected_pt(a, b_end), b_end),
			CgmVec3_sub(CgmCapsule3d_get_projected_pt(b, a_start), a_start),
			CgmVec3_sub(CgmCapsule3d_get_projected_pt(b, a_end), a_end),
		};

		float end_to_projected_pt_distance_squareds[] = {
			CgmVec3_dot(end_to_projected_pt_vecs[0], end_to_projected_pt_vecs[0]),
			CgmVec3_dot(end_to_projected_pt_vecs[1], end_to_projected_pt_vecs[1]),
			CgmVec3_dot(end_to_projected_pt_vecs[2], end_to_projected_pt_vecs[2]),
			CgmVec3_dot(end_to_projected_pt_vecs[3], end_to_projected_pt_vecs[3]),
		};

		for (int i = 0; i < 4; i += 1) {
			if (end_to_projected_pt_distance_squareds[i] < closest_distance_squared) {
				closest_end_to_projected_pt_vec = end_to_projected_pt_vecs[i];
				closest_distance_squared = end_to_projected_pt_distance_squareds[i];
			}
		}
	}

	//
	// all of the separation axis that do not need a cross product to work out.
	CgmVec3 face_axises[] = {
		CgmVec3_norm(a_edge_vector),
		CgmVec3_norm(b_edge_vector),
		CgmVec3_norm(closest_end_to_projected_pt_vec),
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;
	for (uint32_t i = 0; i < 4; i += 1) {
		CgmVec3 axis;
		if (i < 1) {
			axis = edge_axises[i][0];
			CgmVec3 separation_axis = CgmVec3_mul_cross(axis, edge_axises[i][1]);
			float squared_distance = CgmVec3_dot(separation_axis, separation_axis);
			// if the two axis we just cross product together are parallel, we can skip these.
			// this is only because we calculated the closest_end_to_projected_pt_vec which
			// will handle the separation that this cross product was trying to do.
			if (cgm_approx_eq(squared_distance, 0.f)) continue;
			axis = CgmVec3_norm(separation_axis);
		} else {
			axis = face_axises[i - 1];
		}

		//
		// project both of the capsule's start and end onto the axis
		float a_p_start = CgmVec3_dot(a_start, axis);
		float a_p_end = CgmVec3_dot(a_end, axis);
		float a_min = cgm_min(a_p_start, a_p_end) - a->radius;
		float a_max = cgm_max(a_p_start, a_p_end) + a->radius;

		float b_p_start = CgmVec3_dot(b_start, axis);
		float b_p_end = CgmVec3_dot(b_end, axis);
		float b_min = cgm_min(b_p_start, b_p_end) - b->radius;
		float b_max = cgm_max(b_p_start, b_p_end) + b->radius;

		//
		// get smallest distance between the two possible overlapping ends.
		// the smaller is the possible separation or overlap between  the two shapes.
		// the larger one will be the span of the two shapes across the axis.
		float distance = cgm_min(a_max - b_min, b_max - a_min);
		if (distance < 0.f) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		if (distance < min_axis_distance) {
			//
			// we have found a axis with the new minimum separation needed to separate.
			// we need to flip the axis if the capsule B is past capsule A
			min_axis_distance = distance;
			min_axis = a_min + ((a_max - a_min) * 0.5f) < b_min + ((b_max - b_min) * 0.5f) ? CgmVec3_neg(axis) : axis;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_capsule_vs_aabb(CgmCapsule3d* capsule, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_aabb_vs_capsule(aabb, capsule, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_capsule_vs_sphere(CgmCapsule3d* capsule, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_sphere_vs_capsule(sphere, capsule, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_capsule_vs_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	//
	// calculate the start and end positions the capsule.
	// as well as an edge vector that goes from the start to end of the capsule.
	CgmVec3 capsule_offset = CgmVec3_mul_scalar(capsule->direction, capsule->half_length);
	CgmVec3 capsule_start = CgmVec3_add(capsule->center_pos, capsule_offset);
	CgmVec3 capsule_end = CgmVec3_add(capsule->center_pos, CgmVec3_neg(capsule_offset));
	CgmVec3 capsule_edge_vector = CgmVec3_sub(capsule_end, capsule_start);

	//
	// calculate the edge vectors of the triangle.
	// while doing so, find the best axis to test the capsule's start and end
	// with each triangle edge.
	// we do this by first projecting a capsule end point onto an edge. this
	// will give us the closest point on the edge to the capsule end point.
	// we then use this new projected point to find the closest point
	// along the capsule using the same technique.
	// subtracting these two points will give us the separation axis vector we need.
	CgmVec3 triangle_edge_vecs[3];
	CgmVec3 face_axises[7];
	for (uint32_t i = 0; i < 3; i += 1) {
		CgmVec3 start = triangle->arr[i];
		CgmVec3 end = triangle->arr[i + 1 == 3 ? 0 : i + 1];
		CgmVec3 edge = CgmVec3_sub(end, start);
		triangle_edge_vecs[i] = edge;

		float edge_distance_squared = CgmVec3_dot(edge, edge);

		// project capsule_start onto the triangle edge
		float t = CgmVec3_dot(CgmVec3_sub(capsule_start, start), edge) / edge_distance_squared;
		CgmVec3 projected_pt_start = CgmVec3_add(start, CgmVec3_mul_scalar(edge, cgm_clamp(t, 0.f, 1.f)));

		// project capsule_end onto the triangle edge
		t = CgmVec3_dot(CgmVec3_sub(capsule_end, start), edge) / edge_distance_squared;
		CgmVec3 projected_pt_end = CgmVec3_add(start, CgmVec3_mul_scalar(edge, cgm_clamp(t, 0.f, 1.f)));

		// reproject from the points back onto the capsule
		CgmVec3 reprojected_pt_on_capsule_start = CgmCapsule3d_get_projected_pt(capsule, projected_pt_start);
		CgmVec3 reprojected_pt_on_capsule_end = CgmCapsule3d_get_projected_pt(capsule, projected_pt_end);

		// create the separation axis
		face_axises[1 + i * 2] = CgmVec3_norm(CgmVec3_sub(projected_pt_start, reprojected_pt_on_capsule_start));
		face_axises[1 + i * 2 + 1] = CgmVec3_norm(CgmVec3_sub(projected_pt_end, reprojected_pt_on_capsule_end));
	}

	CgmVec3 triangle_face_normal = CgmVec3_norm(CgmVec3_mul_cross(triangle_edge_vecs[0], triangle_edge_vecs[1]));
	face_axises[0] = triangle_face_normal;

	//
	// all combinations of the edges from both shapes paired together.
	// we will use the cross product on these pairs to get the
	// separation axis
	CgmVec3 edge_axises[][2] = {
		{ triangle_edge_vecs[0], capsule_edge_vector },
		{ triangle_edge_vecs[1], capsule_edge_vector },
		{ triangle_edge_vecs[2], capsule_edge_vector },
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;
	for (uint32_t i = 0; i < 9; i += 1) {
		CgmVec3 axis;
		if (i < 3) {
			axis = edge_axises[i][0];
			CgmVec3 separation_axis = CgmVec3_mul_cross(axis, edge_axises[i][1]);
			float squared_distance = CgmVec3_dot(separation_axis, separation_axis);
			// if the two axis we just cross product together are parallel, we can skip these.
			// this is only because we test the projected start and end points of the capsule on each triangle edge.
			if (cgm_approx_eq(squared_distance, 0.f)) continue;
			axis = CgmVec3_norm(separation_axis);
		} else {
			axis = face_axises[i - 3];
		}

		//
		// project the capsule's start and end onto the axis
		float capsule_p_start = CgmVec3_dot(capsule_start, axis);
		float capsule_p_end = CgmVec3_dot(capsule_end, axis);
		float capsule_min = cgm_min(capsule_p_start, capsule_p_end) - capsule->radius;
		float capsule_max = cgm_max(capsule_p_start, capsule_p_end) + capsule->radius;

		// project all 3 vertices of the triangle onto the separation axis
		float triangle_p0 = CgmVec3_dot(triangle->a, axis);
		float triangle_p1 = CgmVec3_dot(triangle->b, axis);
		float triangle_p2 = CgmVec3_dot(triangle->c, axis);
		float triangle_min = cgm_min(cgm_min(triangle_p0, triangle_p1), triangle_p2);
		float triangle_max = cgm_max(cgm_max(triangle_p0, triangle_p1), triangle_p2);

		//
		// get smallest distance between the two possible overlapping ends.
		// the smaller is the possible separation or overlap between the two shapes.
		// the larger one will be the span of the two shapes across the axis.
		float distance = cgm_min(capsule_max - triangle_min, triangle_max - capsule_min);
		if (distance < 0.f) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		if (distance < min_axis_distance) {
			//
			// we have found a axis with the new minimum separation needed to separate.
			// we need to flip the axis if the triangle is in the positive half of the sphere (sphere center is 0 AKA origin).
			min_axis_distance = distance;
			min_axis = capsule_min + ((capsule_max - capsule_min) * 0.5f) < triangle_min + ((triangle_max - triangle_min) * 0.5f) ? CgmVec3_neg(axis) : axis;
		}
	}

	//
	// when using single side triangles, we only want to provide
	// a collision separation when the minimum separation axis
	// is the face normal of the triangle
	if (triangle_single_side) {
		float d = CgmVec3_dot(min_axis, triangle_face_normal);
		if (!cgm_approx_eq(d, 1.f)) {
			return cgm_false;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_capsule_vs_single_side_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_capsule_vs_triangle(capsule, triangle, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_capsule_vs_double_side_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_capsule_vs_triangle(capsule, triangle, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_capsule_vs_mesh(CgmCapsule3d* capsule, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return _cgm_3d_shape_vs_mesh(capsule, mesh, contact_normal_out, contact_distance_out, (CgmBool(*)())cgm_3d_capsule_vs_triangle);
}

CgmBool cgm_3d_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmBool res = cgm_3d_aabb_vs_triangle(aabb, triangle, contact_normal_out, contact_distance_out, triangle_single_side);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmBool res = cgm_3d_sphere_vs_triangle(sphere, triangle, contact_normal_out, contact_distance_out, triangle_single_side);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmBool res = cgm_3d_capsule_vs_triangle(capsule, triangle, contact_normal_out, contact_distance_out, triangle_single_side);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_triangle_vs_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side_a, CgmBool triangle_single_side_b) {
	//
	// compute edge vectors for both of the triangles
	CgmVec3 a_e0 = CgmVec3_sub(a->b, a->a);
	CgmVec3 a_e1 = CgmVec3_sub(a->c, a->b);
	CgmVec3 a_e2 = CgmVec3_sub(a->a, a->c);
	CgmVec3 a_face_normal = CgmVec3_norm(CgmVec3_mul_cross(a_e0, a_e1));

	CgmVec3 b_e0 = CgmVec3_sub(b->b, b->a);
	CgmVec3 b_e1 = CgmVec3_sub(b->c, b->b);
	CgmVec3 b_e2 = CgmVec3_sub(b->a, b->c);
	CgmVec3 b_face_normal = CgmVec3_norm(CgmVec3_mul_cross(b_e0, b_e1));

	//
	// all combinations of the edges from both shapes paired together.
	// we will use the cross product on these pairs to get the
	// separation axis
	CgmVec3 edge_axises[][2] = {
		{ a_e0, b_e0 },
		{ a_e0, b_e1 },
		{ a_e0, b_e2 },

		{ a_e1, b_e0 },
		{ a_e1, b_e1 },
		{ a_e1, b_e2 },

		{ a_e2, b_e0 },
		{ a_e2, b_e1 },
		{ a_e2, b_e2 },
	};

	//
	// all of the separation axis that do not need a cross product to work out.
	CgmVec3 face_axises[] = {
		a_face_normal,
		b_face_normal,
	};

	CgmVec3 min_axis;
	float min_axis_distance = INFINITY;
	for (uint32_t i = 0; i < 11; i += 1) {
		CgmVec3 axis;
		if (i < 9) {
			axis = edge_axises[i][0];
			CgmVec3 separation_axis = CgmVec3_mul_cross(axis, edge_axises[i][1]);
			float squared_distance = CgmVec3_dot(separation_axis, separation_axis);
			// if the two axis we just cross product together are parallel, we need to find
			// a perpendicular axis to test these two edges.
			// we do this by first projecting the triangle A edge start point onto the triangle B edge.
			// this will give us the closest point on the triangle B edge to the triangle A start point.
			// we then use this new projected point to find the closest point
			// back on the triangle A edge using the same technique.
			// subtracting these two points give us for the separation axis vector we need.
			if (cgm_approx_eq(squared_distance, 0.f)) {
				uint32_t a_i = i / 3;
				uint32_t b_i = i % 3;

				// get triangle A edge
				CgmVec3 a_start = a->arr[a_i];
				CgmVec3 a_end = a->arr[a_i + 1 == 3 ? 0 : a_i + 1];

				// get triangle B edge
				CgmVec3 b_start = b->arr[b_i];
				CgmVec3 b_end = b->arr[b_i + 1 == 3 ? 0 : b_i + 1];

				// project A onto B then the result back onto A
				CgmVec3 projected_pt_a_onto_b = cgm_line_get_projected_pt(b_start, b_end, a_start);
				CgmVec3 reprojected_pt_b_onto_a = cgm_line_get_projected_pt(a_start, a_end, projected_pt_a_onto_b);

				separation_axis = CgmVec3_sub(projected_pt_a_onto_b, reprojected_pt_b_onto_a);
			}
			axis = CgmVec3_norm(separation_axis);
		} else {
			axis = face_axises[i - 9];
		}

		// project all 3 vertices of "triangle A" onto the separation axis
		float a_p0 = CgmVec3_dot(a->a, axis);
		float a_p1 = CgmVec3_dot(a->b, axis);
		float a_p2 = CgmVec3_dot(a->c, axis);
		float a_min = cgm_min(a_p0, cgm_min(a_p1, a_p2));
		float a_max = cgm_max(a_p0, cgm_max(a_p1, a_p2));

		// project all 3 vertices of "triangle B" onto the separation axis
		float b_p0 = CgmVec3_dot(b->a, axis);
		float b_p1 = CgmVec3_dot(b->b, axis);
		float b_p2 = CgmVec3_dot(b->c, axis);
		float b_min = cgm_min(b_p0, cgm_min(b_p1, b_p2));
		float b_max = cgm_max(b_p0, cgm_max(b_p1, b_p2));

		// find the smallest distances between the two possible overlaps
		float distance = cgm_min(a_max - b_min, b_max - a_min);
		if (distance < 0.f) {
			//
			// here we have found a separation along this axis
			// so return no collision.
			//
			return cgm_false;
		}

		if (distance > 0.f && distance < min_axis_distance) {
			min_axis_distance = distance;
			min_axis = a_min + ((a_max - a_min) * 0.5f) < b_min + ((b_max - b_min) * 0.5f) ? CgmVec3_neg(axis) : axis;
		}
	}

	//
	// when using single side triangles, we only want to provide
	// a collision separation when the minimum separation axis
	// is the face normal of the triangle
	if (triangle_single_side_a) {
		float d = CgmVec3_dot(min_axis, a_face_normal);
		if (!cgm_approx_eq(d, 1.f)) {
			return cgm_false;
		}
	}
	if (triangle_single_side_b) {
		float d = CgmVec3_dot(min_axis, b_face_normal);
		if (!cgm_approx_eq(d, 1.f)) {
			return cgm_false;
		}
	}

	//
	// we have made it this far so we have a collision.
	if (contact_normal_out) *contact_normal_out = min_axis;
	if (contact_distance_out) *contact_distance_out = min_axis_distance;
	return cgm_true;
}

CgmBool cgm_3d_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmTriangle3d mesh_triangle;
	CgmVec3 best_contact_normal;
	float max_contact_distance = -INFINITY;
	CgmVec3 contact_normal;
	float contact_distance;
	for (uint32_t i = 0; i < mesh->triangles_count; i += 1) {
		mesh_triangle = mesh->triangles[i];
		CgmTriangle3d_offset(&mesh_triangle, mesh->center_pos);
		if (cgm_3d_triangle_vs_triangle(triangle, &mesh_triangle, &contact_normal, &contact_distance, triangle_single_side, cgm_true)) {
			if (!contact_normal_out && !contact_distance_out) return cgm_true;

			if (contact_distance > max_contact_distance) {
				max_contact_distance = contact_distance;
				best_contact_normal = contact_normal;
			}
		}
	}

	if (isinf(max_contact_distance)) return cgm_false;

	if (contact_normal_out) *contact_normal_out = best_contact_normal;
	if (contact_distance_out) *contact_distance_out = contact_distance;
	return cgm_true;
}

CgmBool cgm_3d_single_side_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_aabb(triangle, aabb, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_single_side_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_sphere(triangle, sphere, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_single_side_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_capsule(triangle, capsule, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_single_side_triangle_vs_single_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_triangle(a, b, contact_normal_out, contact_distance_out, cgm_true, cgm_true);
}

CgmBool cgm_3d_single_side_triangle_vs_double_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_triangle(a, b, contact_normal_out, contact_distance_out, cgm_true, cgm_false);
}

CgmBool cgm_3d_single_side_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_mesh(triangle, mesh, contact_normal_out, contact_distance_out, cgm_true);
}


CgmBool cgm_3d_double_side_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_aabb(triangle, aabb, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_double_side_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_sphere(triangle, sphere, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_double_side_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_capsule(triangle, capsule, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_double_side_triangle_vs_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_triangle(a, b, contact_normal_out, contact_distance_out, cgm_false, cgm_true);
}

CgmBool cgm_3d_double_side_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_triangle_vs_mesh(triangle, mesh, contact_normal_out, contact_distance_out, cgm_false);
}

CgmBool cgm_3d_mesh_vs_aabb(CgmMesh* mesh, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_aabb_vs_mesh(aabb, mesh, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_mesh_vs_sphere(CgmMesh* mesh, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_sphere_vs_mesh(sphere, mesh, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_mesh_vs_capsule(CgmMesh* mesh, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out) {
	CgmBool res = cgm_3d_capsule_vs_mesh(capsule, mesh, contact_normal_out, contact_distance_out);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_mesh_vs_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side) {
	CgmBool res = cgm_3d_triangle_vs_mesh(triangle, mesh, contact_normal_out, contact_distance_out, triangle_single_side);
	if (res) {
		if (contact_normal_out) *contact_normal_out = CgmVec3_neg(*contact_normal_out);
	}
	return res;
}

CgmBool cgm_3d_mesh_vs_single_side_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_mesh_vs_triangle(mesh, triangle, contact_normal_out, contact_distance_out, cgm_true);
}

CgmBool cgm_3d_mesh_vs_double_side_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out) {
	return cgm_3d_mesh_vs_triangle(mesh, triangle, contact_normal_out, contact_distance_out, cgm_false);
}

