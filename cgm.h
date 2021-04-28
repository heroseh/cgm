#ifndef CGM_H
#define CGM_H

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef uint8_t CgmBool;
#define cgm_false 0
#define cgm_true 1

#define cgm_assert(cond, message, ...) \
	if (!(cond)) { fprintf(stderr, message"\n", ##__VA_ARGS__); abort(); }

#ifdef CGM_DEBUG_ASSERTIONS
#define cgm_debug_assert cgm_assert
#else
#define cgm_debug_assert(cond, message, ...) (void)(cond)
#endif

#define cgm_debug_uint(v) printf("%s = %u\n", #v, v)
#define cgm_debug_float(v) printf("%s = %f\n", #v, v)
#define cgm_debug_vec2(v) printf("%s = %f, %f\n", #v, (v).x, (v).y)
#define cgm_debug_vec3(v) printf("%s = %f, %f, %f\n", #v, (v).x, (v).y, (v).z)
#define cgm_debug_newline() printf("\n")

typedef union CgmVec2 CgmVec2;
union CgmVec2 {
	struct {
		float x;
		float y;
	};
	float a[2];
};

// ===========================================================================
//
//
// Float Helpers
//
//
// ===========================================================================

#define cgm_epsilon 0.0001

float cgm_min(float a, float b);
float cgm_max(float a, float b);
float cgm_clamp(float v, float min, float max);
float cgm_lerp(float from, float to, float t);
float cgm_lerp_inv(float from, float to, float t);
float cgm_cubic_bezier_curve_interp_1d(CgmVec2 start_anchor, CgmVec2 end_anchor, float ratio);
CgmVec2 cgm_cubic_bezier_curve_interp_2d(CgmVec2 points[4], float ratio);
float cgm_remap(float from_value, float from_min, float from_max, float to_min, float to_max);
CgmBool cgm_approx_eq(float a, float b);

float cgm_sign(float v);
float cgm_round_to_multiple(float v, float multiple);
float cgm_round_up_to_multiple(float v, float multiple);
float cgm_round_down_to_multiple(float v, float multiple);

#define cgm_degrees_to_radians(degrees) ((M_PI / 180.f) * degrees)

typedef struct CgmF16 CgmF16;
struct CgmF16 {
	uint16_t bits;
};

float CgmF16_to_float(CgmF16 v);
CgmF16 CgmF16_from_float(float v);
CgmBool CgmF16_is_nan(CgmF16 v);
CgmBool CgmF16_is_inf(CgmF16 v);

// ===========================================================================
//
//
// Vectors
//
//
// ===========================================================================

typedef union CgmVec2F16 CgmVec2F16;
union CgmVec2F16 {
	struct {
		CgmF16 x;
		CgmF16 y;
	};
	CgmF16 a[2];
};

#define CgmVec2_init(x_, y_) ((CgmVec2){ .x = x_, .y = y_ })
#define CgmVec2_init_even(s) ((CgmVec2){ .x = (s), .y = (s) })
#define CgmVec2_from_f16(v) (CgmVec2){ CgmF16_to_float((v).x), CgmF16_to_float((v).y) }
#define CgmVec2_to_f16(v) (CgmVec2F16){ CgmF16_from_float((v).x), CgmF16_from_float((v).y) }
#define CgmVec2F16_init(x, y) (CgmVec2F16){ CgmF16_from_float(x), CgmF16_from_float(y) }
#define CgmVec2_zero (CgmVec2){0}
#define CgmVec2_inf (CgmVec2){INFINITY, INFINITY}
#define CgmVec2_neg_inf (CgmVec2){-INFINITY, -INFINITY}
CgmBool CgmVec2_eq(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_copysign(CgmVec2 v, CgmVec2 copy);
CgmVec2 CgmVec2_add(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_sub(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_mul(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_div(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_add_scalar(CgmVec2 v, float by);
CgmVec2 CgmVec2_sub_scalar(CgmVec2 v, float by);
CgmVec2 CgmVec2_mul_scalar(CgmVec2 v, float by);
CgmVec2 CgmVec2_div_scalar(CgmVec2 v, float by);

CgmVec2 CgmVec2_neg(CgmVec2 v);
float CgmVec2_len(CgmVec2 v);
CgmVec2 CgmVec2_norm(CgmVec2 v);
float CgmVec2_dot(CgmVec2 a, CgmVec2 b);
float CgmVec2_angle(CgmVec2 v);

CgmVec2 CgmVec2_mul_cross_scalar(CgmVec2 v, float s);
float CgmVec2_mul_cross_vec(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_reflect(CgmVec2 v, CgmVec2 normal);
// absorbs all the force in @param(v) that that goes in the opposite direction as @param(normal)
CgmVec2 CgmVec2_absorb(CgmVec2 v, CgmVec2 normal);
CgmVec2 CgmVec2_rotate(CgmVec2 v, float angle);

#define CgmVec2_up (CgmVec2) { 0.0, -1.0 }
#define CgmVec2_down (CgmVec2) { 0.0, 1.0 }
#define CgmVec2_left (CgmVec2) { -1.0, 0.0 }
#define CgmVec2_right (CgmVec2) { 1.0, 0.0 }
CgmVec2 CgmVec2_perp_left(CgmVec2 v);
CgmVec2 CgmVec2_perp_right(CgmVec2 v);

CgmVec2 CgmVec2_min(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_max(CgmVec2 a, CgmVec2 b);
CgmVec2 CgmVec2_clamp(CgmVec2 v, CgmVec2 min, CgmVec2 max);
CgmVec2 CgmVec2_lerp(CgmVec2 from, CgmVec2 to, CgmVec2 t);
CgmVec2 CgmVec2_sign(CgmVec2 v);
CgmVec2 CgmVec2_round_to_multiple(CgmVec2 v, CgmVec2 multiple);
CgmVec2 CgmVec2_round_up_to_multiple(CgmVec2 v, CgmVec2 multiple);
CgmVec2 CgmVec2_round_down_to_multiple(CgmVec2 v, CgmVec2 multiple);
CgmVec2 CgmVec2_abs(CgmVec2 v);
CgmVec2 CgmVec2_floor(CgmVec2 v);
CgmVec2 CgmVec2_ceil(CgmVec2 v);
CgmVec2 CgmVec2_round(CgmVec2 v);
CgmBool CgmVec2_approx_eq(CgmVec2 a, CgmVec2 b);

typedef union CgmVec3 CgmVec3;
union CgmVec3 {
	struct {
		float x;
		float y;
		float z;
	};
	float a[3];
};

#define CgmVec3_init(x_, y_, z_) ((CgmVec3){ .x = x_, .y = y_, .z = z_ })
#define CgmVec3_init_even(s) ((CgmVec3){ .x = (s), .y = (s), .z = (s) })
#define CgmVec3_zero (CgmVec3){0}
CgmBool CgmVec3_eq(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_copysign(CgmVec3 v, CgmVec3 copy);

CgmVec3 CgmVec3_add(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_sub(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_mul(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_div(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_add_scalar(CgmVec3 v, float by);
CgmVec3 CgmVec3_sub_scalar(CgmVec3 v, float by);
CgmVec3 CgmVec3_mul_scalar(CgmVec3 v, float by);
CgmVec3 CgmVec3_div_scalar(CgmVec3 v, float by);


CgmVec3 CgmVec3_neg(CgmVec3 v);
float CgmVec3_len(CgmVec3 v);
CgmVec3 CgmVec3_norm(CgmVec3 v);
float CgmVec3_dot(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_mul_cross(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_reflect(CgmVec3 v, CgmVec3 normal);
// absorbs all the force in @param(v) that that goes in the opposite direction as @param(normal)
CgmVec3 CgmVec3_absorb(CgmVec3 v, CgmVec3 normal);
CgmVec3 CgmVec3_rotate(CgmVec3 point, CgmVec3 rotation_vector, float angle);
CgmVec3 CgmVec3_clamp_magnitude(CgmVec3 v, float min, float max);

CgmVec3 CgmVec3_perp_left(CgmVec3 v);
CgmVec3 CgmVec3_perp_right(CgmVec3 v);
CgmVec3 CgmVec3_perp_forward(CgmVec3 v);
CgmVec3 CgmVec3_perp_backward(CgmVec3 v);

CgmVec3 CgmVec3_min(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_max(CgmVec3 a, CgmVec3 b);
CgmVec3 CgmVec3_clamp(CgmVec3 v, CgmVec3 min, CgmVec3 max);
CgmVec3 CgmVec3_lerp(CgmVec3 from, CgmVec3 to, CgmVec3 t);
CgmVec3 CgmVec3_sign(CgmVec3 v);
CgmVec3 CgmVec3_round_to_multiple(CgmVec3 v, CgmVec3 multiple);
CgmVec3 CgmVec3_round_up_to_multiple(CgmVec3 v, CgmVec3 multiple);
CgmVec3 CgmVec3_round_down_to_multiple(CgmVec3 v, CgmVec3 multiple);
CgmVec3 CgmVec3_abs(CgmVec3 v);
CgmVec3 CgmVec3_floor(CgmVec3 v);
CgmVec3 CgmVec3_ceil(CgmVec3 v);
CgmVec3 CgmVec3_round(CgmVec3 v);
CgmBool CgmVec3_approx_eq(CgmVec3 a, CgmVec3 b);

typedef union CgmVec4 CgmVec4;
union CgmVec4 {
	struct {
		float x;
		float y;
		float z;
		float w;
	};
	float a[4];
};

#define CgmVec4_init(x_, y_, z_, w_) ((CgmVec4){ .x = x_, .y = y_, .z = z_, .w = w_ })
#define CgmVec4_init_even(s) ((CgmVec4){ .x = (s), .y = (s), .z = (s), .w = (s) })
#define CgmVec4_zero (CgmVec4){0}
CgmBool CgmVec4_eq(CgmVec4 v, CgmVec4 copy);
CgmVec4 CgmVec4_copysign(CgmVec4 v, CgmVec4 copy);

// ===========================================================================
//
//
// Matrices - row major order
//
//
// ===========================================================================

typedef union {
	CgmVec2 row[3];
	float a[6];
} CgmMat3x2;

void CgmMat3x2_identity(CgmMat3x2* out);
void CgmMat3x2_identity_translate(CgmMat3x2* out, CgmVec2 v);
void CgmMat3x2_identity_scale(CgmMat3x2* out, CgmVec2 v);
void CgmMat3x2_identity_rotate(CgmMat3x2* out, float angle);
CgmVec2 CgmMat3x2_row(CgmMat3x2* m, uint32_t row_idx);
CgmVec3 CgmMat3x2_column(CgmMat3x2* m, uint32_t column_idx);

void CgmMat3x2_mul(CgmMat3x2* out, CgmMat3x2* a, CgmMat3x2* b);
CgmVec2 CgmMat3x2_mul_point(CgmMat3x2* m, CgmVec2 pt);
CgmVec2 CgmMat3x2_mul_vector(CgmMat3x2* m, CgmVec2 vector);

typedef union {
	CgmVec4 row[4];
	float a[16];
} CgmMat4x4;

void CgmMat4x4_identity(CgmMat4x4* out);
void CgmMat4x4_identity_scale(CgmMat4x4* out, CgmVec3 v);
void CgmMat4x4_identity_rotate(CgmMat4x4* out, CgmVec3 v, float angle);
void CgmMat4x4_identity_translate(CgmMat4x4* out, CgmVec3 v);

void CgmMat4x4_scale(CgmMat4x4* m, CgmVec3 v);
void CgmMat4x4_rotate(CgmMat4x4* m, CgmVec3 v, float angle);
void CgmMat4x4_translate(CgmMat4x4* m, CgmVec3 v);

void CgmMat4x4_from_3x2(CgmMat4x4* out, CgmMat3x2* m);
CgmVec4 CgmMat4x4_row(CgmMat4x4* m, uint32_t row_idx);
CgmVec4 CgmMat4x4_column(CgmMat4x4* m, uint32_t column_idx);
void CgmMat4x4_ortho(CgmMat4x4* out, float left, float right, float bottom, float top, float near, float far);
void CgmMat4x4_perspective(CgmMat4x4* out, float fovy, float aspect_ratio, float z_near, float z_far);
void CgmMat4x4_mul(CgmMat4x4* out, CgmMat4x4* a, CgmMat4x4* b);
CgmVec3 CgmMat4x4_mul_point(CgmMat4x4* m, CgmVec3 pt);
CgmVec3 CgmMat4x4_mul_vector(CgmMat4x4* m, CgmVec3 v);

// ===========================================================================
//
//
// Collision Shapes 2D
//
//
// ===========================================================================

typedef struct CgmRay2d CgmRay2d;
struct CgmRay2d {
	union {
		struct {
			float x;
			float y;
		};
		CgmVec2 pos;
	};
	CgmVec2 dir;
};

typedef union CgmAabb2d CgmAabb2d;
union CgmAabb2d {
	struct {
		float x;
		float y;
		float ex;
		float ey;
	};
	struct {
		CgmVec2 min;
		CgmVec2 max;
	};
};

#define CgmAabb2d_init(x_, y_, ex_, ey_) (CgmAabb2d){.x = x_, .y = y_, .ex = ex_, .ey = ey_}
#define CgmAabb2d_init_wh(x_, y_, w, h) (CgmAabb2d){.x = x_, .y = y_, .ex = x_ + w, .ey = y_ + h}
CgmVec2 CgmAabb2d_top_right(CgmAabb2d* a);
CgmVec2 CgmAabb2d_bottom_left(CgmAabb2d* a);
float CgmAabb2d_width(CgmAabb2d* a);
float CgmAabb2d_height(CgmAabb2d* a);
CgmVec2 CgmAabb2d_size(CgmAabb2d* a);
CgmVec2 CgmAabb2d_half_size(CgmAabb2d* a);
CgmVec2 CgmAabb2d_center(CgmAabb2d* a);

typedef struct CgmCircle CgmCircle;
struct CgmCircle {
	union {
		struct {
			float x;
			float y;
		};
		CgmVec2 center_pos;
	};
	float radius;
};

typedef struct CgmCapsule2d CgmCapsule2d;
struct CgmCapsule2d {
	CgmVec2 center_pos;
	CgmVec2 direction;
	float half_length;
	float radius;
};

// ===========================================================================
//
//
// Collision Checking 2D
//
//
// ===========================================================================

CgmBool cgm_2d_ray_vs_aabb(CgmRay2d* r, float r_max_distance, CgmAabb2d* a, CgmVec2* contact_normal_out, float* contact_distance_out);
#define cgm_2d_ray_vs_aabb(a, r, r_max_distance, contact_normal_out, contact_distance_out) cgm_2d_aabb_vs_ray(r, r_max_distance, a, contact_normal_out, contact_distance_out)

CgmBool cgm_2d_swept_aabb_vs_aabb(CgmAabb2d* a, CgmVec2 a_dir, float a_distance, CgmAabb2d* b, CgmVec2* contact_normal_out, float* contact_distance_out);
CgmBool cgm_2d_aabb_vs_swept_aabb(CgmAabb2d* a, CgmAabb2d* b, CgmVec2 b_dir, float b_max_distance, CgmVec2* contact_normal_out, float* contact_distance_out);
CgmBool cgm_2d_swept_aabb_vs_swept_aabb(CgmAabb2d* a, CgmVec2 a_dir, float a_distance, CgmAabb2d* b, CgmVec2 b_dir, float b_distance, CgmVec2* contact_normal_out, float* contact_distance_out);
CgmBool cgm_2d_aabb_vs_aabb(CgmAabb2d* a, CgmAabb2d* b, CgmVec2* contact_normal_out, float* contact_distance_out);
CgmBool cgm_2d_aabb_vs_pt(CgmAabb2d* a, CgmVec2 pt, CgmVec2* contact_normal_out, float* contact_distance_out);
#define cgm_2d_pt_vs_aabb(pt, a, contact_normal_out, contact_distance_out) cgm_2d_aabb_vs_pt(a, pt, contact_normal_out, contact_distance_out)

// ===========================================================================
//
//
// Collision Shapes 3D
//
//
// ===========================================================================

typedef struct CgmRay3d CgmRay3d;
struct CgmRay3d {
	CgmVec3 pos;
	CgmVec3 dir;
};
CgmVec3 CgmRay3d_get_projected_pt(CgmRay3d* ray, float ray_max_contact_distance, CgmVec3 pt);

typedef union CgmAabb3d CgmAabb3d;
union CgmAabb3d {
	struct {
		float x;
		float y;
		float z;
		float ex;
		float ey;
		float ez;
	};
	struct {
		CgmVec3 min;
		CgmVec3 max;
	};
};

#define CgmAabb3d_init(x_, y_, z_, ex_, ey_, ez_) ((CgmAabb3d) { .x = (x_), .y = (y_), .z = (z_), .ex = (ex_), .ey = (ey_), .ez = (ez_) })
#define CgmAabb3d_init_min_max(min_, max_) ((CgmAabb3d) { .min = (min_), .max = (max_) })
float CgmAabb3d_width(CgmAabb3d* a);
float CgmAabb3d_height(CgmAabb3d* a);
float CgmAabb3d_depth(CgmAabb3d* a);
CgmVec3 CgmAabb3d_size(CgmAabb3d* a);
CgmVec3 CgmAabb3d_half_size(CgmAabb3d* a);
CgmVec3 CgmAabb3d_center(CgmAabb3d* a);
CgmVec3 CgmAabb3d_clamp_pt(CgmAabb3d* a, CgmVec3 pt);
CgmVec3 CgmAabb3d_closest_pt(CgmAabb3d* a, CgmVec3 pt);
CgmVec3 CgmAabb3d_closest_pt_on_ray(CgmAabb3d* a, CgmVec3 origin, CgmVec3 direction, float min, float max);
CgmVec3 CgmAabb3d_closest_pt_along_ray(CgmAabb3d* a, CgmVec3 origin, CgmVec3 direction, float min, float max);

typedef struct CgmSphere CgmSphere;
struct CgmSphere {
	union {
		struct {
			float x;
			float y;
			float z;
		};
		CgmVec3 center_pos;
	};
	float radius;
};
void CgmSphere_aabb(CgmSphere* sphere, CgmAabb3d* aabb_in_out);

typedef struct CgmCapsule3d CgmCapsule3d;
struct CgmCapsule3d {
	CgmVec3 center_pos;
	CgmVec3 direction;
	float half_length;
	float radius;
};

void CgmCapsule3d_aabb(CgmCapsule3d* capsule, CgmAabb3d* aabb_in_out);
CgmVec3 CgmCapsule3d_get_projected_pt(CgmCapsule3d* capsule, CgmVec3 pt);

typedef union CgmTriangle3d CgmTriangle3d;
union CgmTriangle3d {
	struct {
		CgmVec3 a;
		CgmVec3 b;
		CgmVec3 c;
	};
	CgmVec3 arr[3];
};

void CgmTriangle3d_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb_in_out);
void CgmTriangle3d_offset(CgmTriangle3d* triangle, CgmVec3 by);

CgmVec3 cgm_line_get_projected_pt(CgmVec3 start, CgmVec3 end, CgmVec3 pt);

typedef struct CgmMesh CgmMesh;
struct CgmMesh {
	CgmTriangle3d* triangles;
	CgmVec3 center_pos;
	uint32_t triangles_count;
};

typedef struct CgmHoop3d CgmHoop3d;
struct CgmHoop3d {
	CgmVec3 center_pos;

	// these must be normalized and perpendicular from one another
	CgmVec3 up;
	CgmVec3 right;

	float inner_radius;
	float rim_radius;
};

CgmVec3 CgmHoop3d_closest_pt(CgmHoop3d* hoop, CgmVec3 pt);
void CgmHoop3d_aabb(CgmHoop3d* hoop, CgmAabb3d* aabb_in_out);


// ===========================================================================
//
//
// Collision Checking 3D
//
//
// ===========================================================================

CgmBool cgm_3d_pt_vs_aabb(CgmVec3 pt, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_pt_vs_sphere(CgmVec3 pt, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_pt_vs_capsule(CgmVec3 pt, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);

//
// TODO: test and correct all of these ray functions
//
CgmBool cgm_3d_ray_vs_aabb(CgmRay3d* ray, float ray_max_contact_distance, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_ray_vs_sphere(CgmRay3d* ray, float ray_max_contact_distance, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_ray_vs_capsule(CgmRay3d* ray, float ray_max_contact_distance, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_ray_vs_triangle(CgmRay3d* ray, float ray_max_contact_distance, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_ray_vs_single_side_triangle(CgmRay3d* ray, float ray_max_contact_distance, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_ray_vs_double_side_triangle(CgmRay3d* ray, float ray_max_contact_distance, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_ray_vs_mesh(CgmRay3d* ray, float ray_max_contact_distance, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_aabb_vs_pt(CgmAabb3d* aabb, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_aabb(CgmAabb3d* a, CgmAabb3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_sphere(CgmAabb3d* aabb, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_swept_sphere(CgmAabb3d* aabb, CgmSphere* sphere, CgmVec3 sphere_dir, float sphere_max_contact_distance, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_capsule(CgmAabb3d* aabb, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_aabb_vs_single_side_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_double_side_triangle(CgmAabb3d* aabb, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_aabb_vs_mesh(CgmAabb3d* aabb, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);
// TODO: test and correct AABB vs hoop
CgmBool cgm_3d_aabb_vs_hoop(CgmAabb3d* aabb, CgmHoop3d* hoop, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_sphere_vs_pt(CgmSphere* sphere, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_sphere(CgmSphere* a, CgmSphere* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_aabb(CgmSphere* sphere, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_swept_sphere_vs_aabb(CgmSphere* sphere, CgmVec3 sphere_dir, float sphere_max_contact_distance, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_capsule(CgmSphere* sphere, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_sphere_vs_single_side_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_double_side_triangle(CgmSphere* sphere, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_mesh(CgmSphere* sphere, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_sphere_vs_hoop(CgmSphere* sphere, CgmHoop3d* hoop, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_capsule_vs_pt(CgmCapsule3d* capsule, CgmVec3 pt, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_capsule(CgmCapsule3d* a, CgmCapsule3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_aabb(CgmCapsule3d* capsule, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_sphere(CgmCapsule3d* capsule, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_capsule_vs_single_side_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_double_side_triangle(CgmCapsule3d* capsule, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_mesh(CgmCapsule3d* capsule, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_capsule_vs_hoop(CgmCapsule3d* capsule, CgmHoop3d* hoop, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_triangle_vs_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side_a, CgmBool triangle_single_side_b);
CgmBool cgm_3d_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_triangle_vs_hoop(CgmTriangle3d* triangle, CgmHoop3d* hoop, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);

CgmBool cgm_3d_single_side_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_single_side_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_single_side_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_single_side_triangle_vs_single_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_single_side_triangle_vs_double_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_single_side_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_double_side_triangle_vs_aabb(CgmTriangle3d* triangle, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_double_side_triangle_vs_sphere(CgmTriangle3d* triangle, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_double_side_triangle_vs_capsule(CgmTriangle3d* triangle, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_double_side_triangle_vs_single_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_double_side_triangle_vs_double_side_triangle(CgmTriangle3d* a, CgmTriangle3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_double_side_triangle_vs_mesh(CgmTriangle3d* triangle, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_mesh_vs_aabb(CgmMesh* mesh, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_mesh_vs_sphere(CgmMesh* mesh, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_mesh_vs_capsule(CgmMesh* mesh, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_mesh_vs_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_mesh_vs_hoop(CgmMesh* mesh, CgmHoop3d* hoop, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_mesh_vs_single_side_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_mesh_vs_double_side_triangle(CgmMesh* mesh, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);

CgmBool cgm_3d_hoop_vs_aabb(CgmHoop3d* hoop, CgmAabb3d* aabb, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_sphere(CgmHoop3d* hoop, CgmSphere* sphere, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_capsule(CgmHoop3d* hoop, CgmCapsule3d* capsule, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_triangle(CgmHoop3d* hoop, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out, CgmBool triangle_single_side);
CgmBool cgm_3d_hoop_vs_single_side_triangle(CgmHoop3d* hoop, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_double_side_triangle(CgmHoop3d* hoop, CgmTriangle3d* triangle, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_mesh(CgmHoop3d* hoop, CgmMesh* mesh, CgmVec3* contact_normal_out, float* contact_distance_out);
CgmBool cgm_3d_hoop_vs_hoop(CgmHoop3d* a, CgmHoop3d* b, CgmVec3* contact_normal_out, float* contact_distance_out);

#endif // CGM_H

