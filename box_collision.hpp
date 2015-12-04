#ifndef BOX_COLLISION_H
#define BOX_COLLISION_H
float get_shift(float min_A, float min_B);
bool box_collision(float *min_A, float *max_A, float *min_B, float *max_B);
bool region_collision(float *min_A, float *max_A, float *min_B, float *max_B, float b);
bool adjacent_subvolume(float *min_A, float *max_A, float *min_B, float *max_B);
bool adjacent_region(float *min_A, float *max_A, float *min_B, float *max_B, float b);
#endif /*BOX_COLLISION_H*/
