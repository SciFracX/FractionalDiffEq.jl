window.BENCHMARK_DATA = {
  "lastUpdate": 1712335995857,
  "repoUrl": "https://github.com/SciFracX/FractionalDiffEq.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "52615090+ErikQQY@users.noreply.github.com",
            "name": "Qingyu Qu",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "48be1127311bdb4798db4711a0e422cc68644703",
          "message": "Merge pull request #108 from ErikQQY/qqy/better_bench\n\nAdd better benchmarks",
          "timestamp": "2024-04-06T00:50:44+08:00",
          "tree_id": "92b826bff9fa651e84c7a88f3ce70618e419c553",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/48be1127311bdb4798db4711a0e422cc68644703"
        },
        "date": 1712335993264,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18587804,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 18020466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 18078300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22337024\nallocs=264985\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}