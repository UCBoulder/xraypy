#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>

const double DEG2RAD = 0.0174532925199432957692369076848861271344287188854;


/*
Takes
int: rows
int: cols
float: incident angle degrees
float: pixel size
float: beam center y
float: beam center x
float: detector distance
float: tilt angle degrees
*/


// Function to process a 2D NumPy array and return two 2D NumPy arrays
static PyObject* transform(PyObject* self, PyObject* args) {
    PyArrayObject *data, *weight;
    double incident_angle, pixel_size, beam_center_y, beam_center_x, det_dist, tilt_angle;
    //npy_intp rows, columns;

    if (!PyArg_ParseTuple(args, "O!O!dddddd", &PyArray_Type, &data, &PyArray_Type, &weight,
                          &incident_angle, &pixel_size, &beam_center_y, &beam_center_x, &det_dist, &tilt_angle)) {
        return NULL;
    }

    if (PyArray_NDIM(data) != 2) {
        PyErr_SetString(PyExc_ValueError, "Image array must be a 2 dimensional array.");
        return NULL;
    }

    if (PyArray_NDIM(weight) != 2) {
        PyErr_SetString(PyExc_ValueError, "Input exposure times (weights) must be 2 dimensional array.");
        return NULL;
    }

    npy_intp *shape = PyArray_SHAPE(data);
    npy_intp dim[2] = {shape[0], shape[1]};
    // rows = PyArray_SHAPE(data)[0];
    // columns = PyArray_SHAPE(data)[1];

    if (shape[0] != PyArray_SHAPE(weight)[0] || shape[1] != PyArray_SHAPE(weight)[1]) {
        PyErr_SetString(PyExc_ValueError, "The data array and the exposure times must be the same sized arrays");
        return NULL;
    }
    
    // Convert angles from degrees to radians
    incident_angle *= DEG2RAD;
    tilt_angle *= DEG2RAD;

    // convert poni file parameters to pixel space
    beam_center_x /= pixel_size;
    beam_center_y /= -pixel_size;
    beam_center_y += (double)shape[0];
    det_dist /= pixel_size;  

    // Create x array
    npy_intp x_dim[2] = {1, shape[1]};
    PyObject* x_obj = PyArray_EMPTY(2, x_dim, NPY_DOUBLE, 0);
    if (x_obj == NULL) {
        return NULL;
    }
    int* x_data = (int*)PyArray_DATA((PyArrayObject*)x_obj);
    for (size_t ii = 0; ii < shape[1]; ++ii) {
        x_data[ii] = beam_center_x - (double)ii;
    }

    // Create y array
    npy_intp y_dim[2] = {shape[1]};
    PyObject* y_obj = PyArray_EMPTY(2, y_dim, NPY_DOUBLE, 0);
    if (y_obj == NULL) {
        Py_DECREF(x_obj);
        return NULL;
    }
    int* y_data = (int*)PyArray_DATA((PyArrayObject*)y_obj);
    for (size_t ii = 0; ii < shape[0]; ++ii) {
        y_data[ii] = beam_center_y - (double)ii;
    }

    if (tilt_angle != 0.0) {
        double cos_tilt = std::cos(tilt_angle);
        double sin_tilt = std::sin(tilt_angle);

        //y_rot = y_lab * tilt_cos - z_lab * tilt_sin
        //z_rot = z_lab * tilt_cos + y_lab * tilt_sin
    }







    // Perform element-wise multiplication using np.multiply
    PyObject* z_obj = PyArray_EMPTY(2, dim, NPY_INT, 0);
    if (z_obj == NULL) {
        Py_DECREF(x_obj);
        Py_DECREF(y_obj);
        return NULL;
    }
    PyObject* multiply_result = PyObject_CallFunctionObjArgs((PyObject*)&PyArray_API[56], x_obj, y_obj, NULL); // np.multiply
    if (multiply_result == NULL) {
        Py_DECREF(x_obj);
        Py_DECREF(y_obj);
        Py_DECREF(z_obj);
        return NULL;
    }
    PyArrayObject* z_array = (PyArrayObject*)z_obj;
    memcpy(PyArray_DATA(z_array), PyArray_DATA((PyArrayObject*)multiply_result), PyArray_NBYTES(z_array));

    // Create output arrays
    npy_intp dims[3] = {2, shape[0], shape[1]};
    PyObject *r = PyArray_EMPTY(3, dims, NPY_DOUBLE, 0);

    // Clean up
    Py_DECREF(x_obj);
    Py_DECREF(y_obj);
    Py_DECREF(multiply_result);




















    // Processing logic here...
    
    // Return the output arrays
    return Py_BuildValue("OO", output_array1, output_array2);
}

// Method definitions
static PyMethodDef methods[] = {
    {"transform", transform, METH_VARARGS, "Process a TIFF image of a thin film exposure."},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "mymodule",
    NULL,
    -1,
    methods
};

// Module initialization function
PyMODINIT_FUNC PyInit_mymodule(void) {
    import_array();  // Initialize NumPy
    return PyModule_Create(&module);
}

