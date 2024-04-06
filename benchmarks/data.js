window.BENCHMARK_DATA = {
  "lastUpdate": 1712381581283,
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
      },
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
          "id": "37963e0f168edc14574c98be212fcc83bdd31355",
          "message": "Merge pull request #109 from ErikQQY/qqy/triger_bench\n\nTrigger benchmarks",
          "timestamp": "2024-04-06T12:58:39+08:00",
          "tree_id": "808e5a17757f6632791169aa8d1c15bc9683f7c7",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/37963e0f168edc14574c98be212fcc83bdd31355"
        },
        "date": 1712379651197,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17310408.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16730695.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16781819.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22337024\nallocs=264985\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "302ce2a3de9d5c0f72d40c2ee1d4299bc55c1a17",
          "message": "Merge pull request #110 from ErikQQY/qqy/better_ben\n\nTest benchmarks",
          "timestamp": "2024-04-06T13:11:05+08:00",
          "tree_id": "5e438e322d22539c126c00d7af0c32a3bb3a264b",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/302ce2a3de9d5c0f72d40c2ee1d4299bc55c1a17"
        },
        "date": 1712380402650,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18292251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17864163.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419536\nallocs=264997\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 17813185.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "a101dbf65fb6651424cec395f5861ed45d61fbed",
          "message": "Merge pull request #111 from ErikQQY/qqy/test_bench\n\nTest CI benchmarks",
          "timestamp": "2024-04-06T13:30:43+08:00",
          "tree_id": "0967cfb08441a2e61a4a0cb00faa22b8423f7403",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/a101dbf65fb6651424cec395f5861ed45d61fbed"
        },
        "date": 1712381578766,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17007258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16916692,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16571171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}