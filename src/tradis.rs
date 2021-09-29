pub struct Trajectory {
    pub(crate) coords: Vec<Coord3D>,
    pub(crate) idx: usize,
}

impl std::ops::Index<usize> for Trajectory {
    type Output = Coord3D;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.coords[idx]
    }
}

impl std::ops::IndexMut<usize> for Trajectory {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.coords[idx]
    }
}

impl Trajectory {
    fn trim(&mut self, start: f64, end: f64) {
        let trj = self.coords.to_owned();
        let mut start_idx: usize = 0;
        let mut end_idx: usize = trj.len() - 1;
        let mut new_start: Coord3D = trj[start_idx].to_owned();
        let mut new_end: Coord3D = trj[end_idx].to_owned();
        let mut p: Coord3D;
        let mut q: Coord3D;
        for i in 0..(trj.len() - 1) {
            p = trj[i].to_owned();
            q = trj[i + 1].to_owned();
            if (p.t < start) & (start < q.t) {
                new_start = interpolate(start, &p, &q);
                start_idx = i;
            };
            if (p.t < end) & (end < q.t) {
                new_end = interpolate(end, &p, &q);
                end_idx = i + 1;
            };
        }

        let mut trimmed: Vec<Coord3D> = trj[start_idx..(end_idx + 1)].to_vec();
        trimmed[0] = new_start;
        trimmed[end_idx - start_idx] = new_end;
        self.coords = trimmed;
    }
    fn len(&self) -> usize {
        self.coords.len()
    }
    fn next_line(&mut self) -> Line3D {
        let trj = &self.coords;
        let idx = self.idx;
        let nxt_idx = if idx < trj.len() - 2 { idx + 1 } else { idx };
        Line3D {
            start: trj[idx].to_owned(),
            end: trj[idx + 1].to_owned(),
            trj_idx: nxt_idx,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Coord3D {
    pub x: f64,
    pub y: f64,
    pub t: f64,
}

struct Line3D {
    start: Coord3D,
    end: Coord3D,
    trj_idx: usize,
}

#[derive(Debug)]
struct Coord2D {
    x: f64,
    y: f64,
}

struct Line2D {
    start: Coord2D,
    end: Coord2D,
}

trait Project {
    fn get_projection(&self, axis: usize) -> Line2D;
}

impl Project for Line3D {
    fn get_projection(&self, axis: usize) -> Line2D {
        let (start_x, end_x) = if axis == 0 {
            (self.start.x, self.end.x)
        } else {
            (self.start.y, self.end.y)
        };

        Line2D {
            start: Coord2D {
                x: start_x,
                y: self.start.t,
            },
            end: Coord2D {
                x: end_x,
                y: self.end.t,
            },
        }
    }
}

pub fn tradis(mut trj_a: Trajectory, mut trj_b: Trajectory) -> f64 {
    let span = get_common_time_span(&trj_a, &trj_b);
    if span == None {
        return f64::INFINITY;
    }
    let (start, end) = span.unwrap();
    trj_a.trim(start, end);
    trj_b.trim(start, end);
    let mut line3d_a: Line3D = trj_a.next_line();
    let mut line3d_b: Line3D = trj_b.next_line();
    let mut areas: [f64; 2] = [0.0, 0.0];
    for (axis, area) in areas.iter_mut().enumerate() {
        while (line3d_a.trj_idx < trj_a.len() - 3) & (line3d_b.trj_idx < trj_b.len() - 3) {
            let line2d_a: Line2D = line3d_a.get_projection(axis);
            let line2d_b: Line2D = line3d_b.get_projection(axis);
            let p: Option<Coord2D> = intersection_point(&line2d_a, &line2d_b);
            match p {
                Some(r) => {
                    *area += calc_area(&line2d_a.start, &line2d_b.start, &r)
                        + calc_area(&r, &line2d_a.end, &line2d_b.end);
                    line3d_a = trj_a.next_line();
                    line3d_b = trj_b.next_line();
                }
                None => {
                    if line2d_a.end.y < line2d_b.end.y {
                        *area += calc_area(&line2d_a.start, &line2d_b.start, &line2d_a.end);
                        line3d_a = trj_a.next_line();
                    } else {
                        *area += calc_area(&line2d_a.start, &line2d_b.start, &line2d_b.end);
                        line3d_b = trj_b.next_line();
                    };
                }
            }
        }
        let line2d_a: Line2D = line3d_a.get_projection(axis);
        let line2d_b: Line2D = line3d_b.get_projection(axis);
        let p: Option<Coord2D> = intersection_point(&line2d_a, &line2d_b);
        match p {
            Some(r) => {
                *area += calc_area(&line2d_a.start, &line2d_b.start, &r)
                    + calc_area(&r, &line2d_a.end, &line2d_b.end);
            }
            None => {
                if line2d_a.end.y < line2d_b.end.y {
                    *area += calc_area(&line2d_a.start, &line2d_b.start, &line2d_a.end);
                    line3d_a = trj_a.next_line();
                } else {
                    *area += calc_area(&line2d_a.start, &line2d_b.start, &line2d_b.end);
                    line3d_b = trj_b.next_line();
                };
            }
        }
    }
    areas[0] + areas[1]
}

fn calc_area(p: &Coord2D, q: &Coord2D, r: &Coord2D) -> f64 {
    0.5 * (p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y)).abs()
}

fn get_common_time_span(trj_a: &Trajectory, trj_b: &Trajectory) -> Option<(f64, f64)> {
    let a_len = trj_a.len() - 1;
    let b_len = trj_b.len() - 1;
    let start: f64 = if trj_a[0].t > trj_b[0].t {
        trj_a[0].t
    } else {
        trj_b[0].t
    };
    let end: f64 = if trj_a[a_len].t < trj_b[b_len].t {
        trj_a[a_len].t
    } else {
        trj_b[b_len].t
    };
    if end - start < 0.0 {
        None
    } else {
        Some((start, end))
    }
}

fn interpolate(t: f64, p: &Coord3D, q: &Coord3D) -> Coord3D {
    let dx = (q.x - p.x) / (q.t - p.t);
    let dy = (q.y - p.y) / (q.t - p.t);
    Coord3D {
        x: p.x + (q.t - t) * dx,
        y: p.y + (q.t - t) * dy,
        t,
    }
}

fn intersection_point(p: &Line2D, q: &Line2D) -> Option<Coord2D> {
    let (x1, y1, x2, y2) = (p.start.x, p.start.y, p.end.x, p.end.y);
    let (x3, y3, x4, y4) = (q.start.x, q.start.y, q.end.x, q.end.y);
    let denom: f64 = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    let t: f64 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
    let u: f64 = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom;
    if (0.0..=1.0).contains(&t) {
        return Some(Coord2D {
            x: x1 + t * (x2 - x1),
            y: y1 + t * (y2 - y1),
        });
    } else if (0.0..=1.0).contains(&u) {
        return Some(Coord2D {
            x: x3 + u * (x4 - x3),
            y: y3 + u * (y4 - y3),
        });
    };
    None
}

#[cfg(test)]
mod tradis_test {
    use super::*;

    #[test]
    fn intersect_simple() {
        let p1 = Coord2D { x: 0.0, y: 0.0 };
        let p2 = Coord2D { x: 2.0, y: 2.0 };
        let q1 = Coord2D { x: 2.0, y: 0.0 };
        let q2 = Coord2D { x: 0.0, y: 2.0 };
        let r = intersection_point(
            &Line2D { start: p1, end: p2 },
            &Line2D { start: q1, end: q2 },
        );
        let r = r.unwrap();
        assert_eq!((r.x, r.y), (1.0, 1.0))
    }

    #[test]
    fn intersect_on_line() {
        let p1 = Coord2D { x: 0.0, y: 0.0 };
        let p2 = Coord2D { x: 2.0, y: 2.0 };
        let q1 = Coord2D { x: 2.0, y: 0.0 };
        let q2 = Coord2D { x: 1.0, y: 1.0 };
        let r = intersection_point(
            &Line2D { start: p1, end: p2 },
            &Line2D { start: q1, end: q2 },
        );
        let r = r.unwrap();
        assert_eq!((r.x, r.y), (1.0, 1.0))
    }

    #[test]
    fn intersect_parallel() {
        let p1 = Coord2D { x: 0.0, y: 0.0 };
        let p2 = Coord2D { x: 0.0, y: 2.0 };
        let q1 = Coord2D { x: 2.0, y: 0.0 };
        let q2 = Coord2D { x: 2.0, y: 2.0 };
        let r = intersection_point(
            &Line2D { start: p1, end: p2 },
            &Line2D { start: q1, end: q2 },
        );
        assert!(Option::is_none(&r))
    }

    #[test]
    fn simple_interpolation() {
        let p: Coord3D = Coord3D {
            x: 0.0,
            y: 0.0,
            t: 0.0,
        };
        let q: Coord3D = Coord3D {
            x: 2.0,
            y: 2.0,
            t: 2.0,
        };
        let ans: Coord3D = Coord3D {
            x: 1.0,
            y: 1.0,
            t: 1.0,
        };
        let res = interpolate(1.0, &p, &q);
        assert_eq!((res.x, res.y, res.t), (ans.x, ans.y, ans.t));
    }
    #[test]
    fn simple_interpolation2() {
        let p: Coord3D = Coord3D {
            x: 0.0,
            y: 0.0,
            t: 0.0,
        };
        let q: Coord3D = Coord3D {
            x: 1.0,
            y: 1.0,
            t: 1.0,
        };
        let ans: Coord3D = Coord3D {
            x: 0.5,
            y: 0.5,
            t: 0.5,
        };
        let res = interpolate(0.5, &q, &p);
        assert_eq!((res.x, res.y, res.t), (ans.x, ans.y, ans.t));
    }
    #[test]
    fn offset_interpolation() {
        let p: Coord3D = Coord3D {
            x: 2.0,
            y: 2.0,
            t: 2.0,
        };
        let q: Coord3D = Coord3D {
            x: 4.0,
            y: 4.0,
            t: 4.0,
        };
        let ans: Coord3D = Coord3D {
            x: 3.0,
            y: 3.0,
            t: 3.0,
        };
        let res = interpolate(3.0, &p, &q);
        assert_eq!((res.x, res.y, res.t), (ans.x, ans.y, ans.t));
    }
    #[test]
    fn single_dim_interpolation() {
        let p: Coord3D = Coord3D {
            x: 2.0,
            y: 2.0,
            t: 2.0,
        };
        let q: Coord3D = Coord3D {
            x: 4.0,
            y: 2.0,
            t: 4.0,
        };
        let ans: Coord3D = Coord3D {
            x: 3.0,
            y: 2.0,
            t: 3.0,
        };
        let res = interpolate(3.0, &p, &q);
        assert_eq!((res.x, res.y, res.t), (ans.x, ans.y, ans.t));
    }
    #[test]
    fn common_time_span() {
        let trj_a = Trajectory {
            coords: vec![
                Coord3D {
                    x: 0.0,
                    y: 0.0,
                    t: 0.0,
                },
                Coord3D {
                    x: 1.0,
                    y: 1.0,
                    t: 1.0,
                },
                Coord3D {
                    x: 2.0,
                    y: 2.0,
                    t: 2.0,
                },
            ],
            idx: 0,
        };
        let trj_b = Trajectory {
            coords: vec![
                Coord3D {
                    x: 0.0,
                    y: 0.0,
                    t: 0.5,
                },
                Coord3D {
                    x: 1.0,
                    y: 1.0,
                    t: 1.5,
                },
                Coord3D {
                    x: 2.0,
                    y: 2.0,
                    t: 2.5,
                },
            ],
            idx: 0,
        };
        let (start, end) = get_common_time_span(&trj_a, &trj_b).unwrap();
        assert!(start < end, "Start is after end");
        assert!((start - end) < 0.00001, "Start equals end");
        assert!((start - 0.5).abs() < 0.0001, "Wrong start time");
        assert!((end - 2.0).abs() < 0.0001, "Wrong end time");
    }
    #[test]
    fn common_time_span2() {
        let trj_a = Trajectory {
            coords: vec![
                Coord3D {
                    x: 0.0,
                    y: 0.0,
                    t: 0.0,
                },
                Coord3D {
                    x: 0.1,
                    y: 0.1,
                    t: 0.1,
                },
                Coord3D {
                    x: 0.2,
                    y: 0.2,
                    t: 0.2,
                },
                Coord3D {
                    x: 1.0,
                    y: 1.0,
                    t: 1.0,
                },
                Coord3D {
                    x: 1.2,
                    y: 1.2,
                    t: 1.2,
                },
                Coord3D {
                    x: 2.0,
                    y: 2.0,
                    t: 2.0,
                },
            ],
            idx: 0,
        };
        let trj_b = Trajectory {
            coords: vec![
                Coord3D {
                    x: 0.0,
                    y: 0.0,
                    t: 0.5,
                },
                Coord3D {
                    x: 1.0,
                    y: 1.0,
                    t: 1.5,
                },
                Coord3D {
                    x: 2.0,
                    y: 2.0,
                    t: 2.5,
                },
            ],
            idx: 0,
        };
        let (start, end) = get_common_time_span(&trj_a, &trj_b).unwrap();
        assert!(start < end, "Start is after end");
        assert!((start - 0.5).abs() < 0.0001);
        assert!((end - 2.0).abs() < 0.0001);
    }

    #[test]
    fn area_simple() {
        let p = Coord2D { x: 0.0, y: 0.0 };
        let q = Coord2D { x: 1.0, y: 0.0 };
        let r = Coord2D { x: 0.0, y: 1.0 };
        assert!((calc_area(&p, &q, &r) - 0.5).abs() < 0.0001)
    }

    #[test]
    fn area_simple2() {
        let p = Coord2D { x: 1.0, y: 1.0 };
        let q = Coord2D { x: 2.0, y: 1.0 };
        let r = Coord2D { x: 1.0, y: 2.0 };
        assert!((calc_area(&p, &q, &r) - 0.5).abs() < 0.0001)
    }
}
