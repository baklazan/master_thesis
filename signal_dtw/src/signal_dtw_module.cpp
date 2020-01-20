#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <local_alignment.h>


namespace py = pybind11;

PYBIND11_MODULE(signal_dtw, m) {
  m.doc() = "A C++ implementation of signal dtw";

  m.def("local_alignment", &local_alignment, "",
        py::arg("signal"), py::arg("chiron_predictions"), py::arg("min_distance"),
        py::arg("max_distance"), py::arg("score_for_moving"), 
        py::arg("max_speed_ratio"), py::arg("event_threshold"),
        py::arg("min_lookahead"), py::arg("log_filename")="");
    

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
