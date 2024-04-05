cimport numpy as np
from libcpp.vector cimport vector
from libc.math cimport sin, cos, asin, sqrt


# cdef extern from "numpy/arrayobject.h":
#     np.ndarray[double, ndim=2] PyArray_EMPTY(int nd, np.npy_intp* dims, int typenum, int fortran, np.npy_intp* strides)



cdef float sign(double x):
    return float((x > 0) - (x < 0))


cdef class Transform:
    cdef double incident_angle
    cdef double tilt_angle
    cdef int[2] shape
    cdef double beam_px_y
    cdef double beam_px_x
    cdef double det_dist_px
    cdef np.ndarray[double, ndim=3] r_xy
    cdef np.ndarray[double, ndim=2] r_z

    def __cinit__(self, float incident_angle_degrees, int rows, int columns,
                  float pixel_size, float beam_center_y, float beam_center_x,
                  float det_dist, float tilt_angle_degrees=0):
        self.incident_angle = incident_angle_degrees * np.pi / 180.
        if tilt_angle_degrees:
            self.tilt_angle = tilt_angle_degrees * np.pi / 180.
        else:
            self.tilt_angle = 0.
        self.shape[0] = rows
        self.shape[1] = columns
        self.beam_px_x = beam_center_x / pixel_size
        self.beam_px_y = <float>rows - beam_center_y / pixel_size
        self.det_dist_px = det_dist / pixel_size
        self.r = np.PyArray_EMPTY(3, (2, rows, columns), np.NPY_DOUBLE, 0)
        self.calc_transformed_r()

    cdef void calc_transformed_r(self):
        cdef double y_lab
        cdef double z_lab
        cdef double tilt_cos
        cdef double tilt_sin
        cdef double sum_dsq_ysq
        cdef double alpha_scattered
        cdef double sin_phi_scattered
        cdef double y_rot
        cdef double z_rot
        cdef double cos_alpha
        cdef double cos_phi
        cdef double q_xy_sq
        cdef double q_xy
        cdef double q_z
        cdef double q_z_sq
        cdef double q_sq
        cdef double q_scaler

        cdef double cos_incid = cos(self.incident_angle)

        for x in range(self.shape[1]):
            for y in range(self.shape[0]):
                y_lab = self.beam_px_x - x
                z_lab = self.beam_px_y - y
                if self.tilt_angle:
                    tilt_cos = cos(self.tilt_angle)
                    tilt_sin = sin(self.tilt_angle)
                    y_rot = y_lab * tilt_cos - z_lab * tilt_sin
                    z_rot = z_lab * tilt_cos + y_lab * tilt_sin
                    sum_dsq_ysq = self.det_dist_px * self.det_dist_px + y_rot * y_rot
                    alpha_scattered = asin(z_rot / sqrt(sum_dsq_ysq + z_rot * z_rot)) - self.incident_angle
                    sin_phi_scattered = y_rot / sqrt(sum_dsq_ysq)
                else:
                    sum_dsq_ysq = self.det_dist_px * self.det_dist_px + y_lab * y_lab
                    alpha_scattered = asin(z_lab / sqrt(sum_dsq_ysq + z_lab * z_lab)) - self.incident_angle
                    sin_phi_scattered = y_lab / sqrt(sum_dsq_ysq)
                cos_alpha = cos(alpha_scattered)
                cos_phi = sqrt(1 - sin_phi_scattered * sin_phi_scattered)
                q_xy_sq = cos_alpha * cos_alpha + cos_incid * cos_incid - (2.0 * cos_incid) * cos_alpha * cos_phi
                q_xy = sqrt(q_xy_sq) * sign(y_lab)
                q_z = sin(alpha_scattered) + sin(self.incident_angle)
                q_z_sq = q_z * q_z
                q_sq = q_xy_sq + q_z_sq
                q_scaler = self.det_dist_px * sqrt(0.5 + 1.0 / (2.0 - q_sq))
                self.r_xy[] = q_scaler * q_xy
                self.r_z[] = q_scaler * q_z



    cpdef (np.ndarray[np.float64_t, ndim=2], np.ndarray[np.float64_t, ndim=2], float, float) transform_image(self, np.ndarray[np.float64_t, ndim=2] data, np.ndarray[np.float64_t, ndim=2] weight):
        cdef float beam_center_xy = np.max(self.r_xy)
        cdef float beam_center_z = np.max(self.r_z)
        cdef int transformed_rows = <int>(np.ceil(beam_center_z - np.min(self.r_z))) + 1
        cdef int transformed_columns = <int>(np.ceil(beam_center_xy - np.min(self.r_xy))) + 1

        cdef np.ndarray[np.float64_t, ndim=2] xy_det_origin = beam_center_xy - self.r_xy
        cdef np.ndarray[np.float64_t, ndim=2] z_det_origin = beam_center_z - self.r_z
        cdef np.ndarray[np.float64_t, ndim=2] xy_det_floor = np.floor(xy_det_origin)
        cdef np.ndarray[np.float64_t, ndim=2] z_det_floor = np.floor(z_det_origin)
        cdef np.ndarray[int, ndim=2] x_pixel_loc = xy_det_floor.astype(int)
        cdef np.ndarray[int, ndim=2] y_pixel_loc = z_det_floor.astype(int)
        cdef np.ndarray[np.float64_t, ndim=2] remainder_xy = xy_det_origin - xy_det_floor
        cdef np.ndarray[np.float64_t, ndim=2] remainder_z = z_det_origin - z_det_floor
        cdef np.ndarray[np.float64_t, ndim=2] remainder_xy_comp = 1.0 - remainder_xy
        cdef np.ndarray[np.float64_t, ndim=2] remainder_z_comp = 1.0 - remainder_z
        cdef np.ndarray[np.float64_t, ndim=2] current_pixel_weight = remainder_xy_comp * remainder_z_comp
        cdef np.ndarray[np.float64_t, ndim=2] x_neighbor_weight = remainder_xy * remainder_z_comp
        cdef np.ndarray[np.float64_t, ndim=2] y_neighbor_weight = remainder_xy_comp * remainder_z
        cdef np.ndarray[np.float64_t, ndim=2] diag_neighbor_weight = remainder_xy * remainder_z

        cdef np.ndarray[np.float64_t, ndim=2] data_t = np.zeros((transformed_rows, transformed_columns), dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=2] weights_t = np.zeros((transformed_rows, transformed_columns), dtype=np.float64)
        for rr in range(self.shape[0]):
            for cc in range(self.shape[1]):
                col = x_pixel_loc[rr, cc]
                row = y_pixel_loc[rr, cc]
                data_t[row, col] += data[rr, cc] * current_pixel_weight[rr, cc]
                weights_t[row, col] += weight[rr, cc] * current_pixel_weight[rr, cc]
                data_t[row, col + 1] += data[rr, cc] * x_neighbor_weight[rr, cc]
                weights_t[row, col + 1] += weight[rr, cc] * x_neighbor_weight[rr, cc]

                data_t[row + 1, col] += data[rr, cc] * y_neighbor_weight[rr, cc]
                weights_t[row + 1, col] += weight[rr, cc] * y_neighbor_weight[rr, cc]

                data_t[row + 1, col + 1] += data[rr, cc] * diag_neighbor_weight[rr, cc]
                weights_t[row + 1, col + 1] += weight[rr, cc] * diag_neighbor_weight[rr, cc]
        return data_t, weights_t, beam_center_z, beam_center_xy

