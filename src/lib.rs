use numpy::PyReadonlyArray2;
use pyo3::{
    prelude::{pymodule, PyModule, PyResult, Python},
    types::PyFloat,
};
use tradis::{Coord3D, Trajectory};
mod tradis;

#[pymodule]
#[pyo3(name = "tradis")]
fn rust_ext(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    fn tradispy(a: tradis::Trajectory, b: tradis::Trajectory) -> f64 {
        tradis::tradis(a, b)
    }

    // wrapper of `axpy`
    #[pyfn(m)]
    #[pyo3(name = "tradis")]
    fn tradispy_py<'py>(
        py: Python<'py>,
        t1: PyReadonlyArray2<f64>,
        t2: PyReadonlyArray2<f64>,
    ) -> &'py PyFloat {
        let a = t1.as_array();
        let b = t2.as_array();
        assert_eq!(a.shape()[1], 3, "Shape of array must be (_,3)");
        assert_eq!(b.shape()[1], 3, "Shape of array must be (_,3)");
        let res: f64 = tradispy(to_trajectory(&t1), to_trajectory(&t2));
        PyFloat::new(py, res)
    }
    Ok(())
}

fn to_trajectory(trj: &PyReadonlyArray2<f64>) -> Trajectory {
    let coords = trj
        .as_array()
        .outer_iter()
        .map(|c| Coord3D {
            x: c[0].to_owned(),
            y: c[1].to_owned(),
            t: c[2].to_owned(),
        })
        .collect();
    Trajectory { coords, idx: 0 }
}
