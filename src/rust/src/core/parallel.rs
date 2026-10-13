//! Deterministic parallel map: item `i` is computed by `f(i)` and stored in
//! slot `i`, whatever the number of threads.

pub fn par_map<T, F>(n_items: usize, n_threads: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync,
{
    let n_threads = n_threads.clamp(1, n_items.max(1));
    if n_threads == 1 {
        return (0..n_items).map(f).collect();
    }
    let per = (n_items + n_threads - 1) / n_threads;
    let f = &f;
    let chunks: Vec<Vec<T>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n_threads)
            .map(|w| {
                let lo = (w * per).min(n_items);
                let hi = ((w + 1) * per).min(n_items);
                scope.spawn(move || (lo..hi).map(f).collect::<Vec<T>>())
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("worker thread panicked")).collect()
    });
    chunks.into_iter().flatten().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_is_kept() {
        for threads in 1..5 {
            assert_eq!(par_map(11, threads, |i| i * i), (0..11).map(|i| i * i).collect::<Vec<_>>());
        }
        assert!(par_map(0, 3, |i| i).is_empty());
    }
}
