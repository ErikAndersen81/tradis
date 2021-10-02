pub mod util;

use numpy::PyReadonlyArray2;
use pyo3::{
    prelude::{pymodule, PyModule, PyResult, Python},
    types::PyFloat,
};

pub use util::{distance, Trajectory};

#[pymodule]
#[pyo3(name = "tradis")]
fn rust_ext(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    fn tradispy(a: Vec<[f64; 3]>, b: Vec<[f64; 3]>) -> f64 {
        let trj_a = Trajectory::from_array(a);
        let trj_b = Trajectory::from_array(b);
        let dist = distance(&trj_a, &trj_b);
        match dist {
            Some(dist) => dist,
            None => f64::NAN,
        }
    }

    // wrapper of `axpy`
    #[pyfn(m)]
    #[pyo3(name = "tradis")]
    fn tradispy_py<'py>(
        py: Python<'py>,
        trj_a: PyReadonlyArray2<f64>,
        trj_b: PyReadonlyArray2<f64>,
    ) -> &'py PyFloat {
        let trj_a = trj_a.as_array();
        let trj_b = trj_b.as_array();
        assert_eq!(trj_a.shape()[1], 3, "Shape of array must be (_,3)");
        assert_eq!(trj_b.shape()[1], 3, "Shape of array must be (_,3)");
        assert!(
            trj_a.len() > 1,
            "Trajectories must have more than one point"
        );
        assert!(
            trj_b.len() > 1,
            "Trajectories must have more than one point"
        );

        let trj_a = trj_a
            .outer_iter()
            .map(|c| [c[0], c[1], c[2]])
            .collect::<Vec<[f64; 3]>>();
        let trj_b = trj_b
            .outer_iter()
            .map(|c| [c[0], c[1], c[2]])
            .collect::<Vec<[f64; 3]>>();
        let res: f64 = tradispy(trj_a, trj_b);
        PyFloat::new(py, res)
    }
    Ok(())
}
