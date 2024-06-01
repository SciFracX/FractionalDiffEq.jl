window.BENCHMARK_DATA = {
  "lastUpdate": 1717249147180,
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
          "id": "61a089ca536331e7dfa346e1d3554aa8474da65d",
          "message": "Merge pull request #112 from ErikQQY/qqy/test_bench\n\nTrigger again",
          "timestamp": "2024-06-01T20:39:28+08:00",
          "tree_id": "ebf3007828d12bd460e69e68a28e29deab339737",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/61a089ca536331e7dfa346e1d3554aa8474da65d"
        },
        "date": 1717245718580,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17357485,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22961808\nallocs=265810\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16803226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22417184\nallocs=264853\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16816974.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22334672\nallocs=264841\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "4c10b03f697e027bf99839b130a372a730392d2e",
          "message": "Add formatter\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T20:41:21+08:00",
          "tree_id": "89bf8adc6516b1ca4ee251500427894761c01381",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/4c10b03f697e027bf99839b130a372a730392d2e"
        },
        "date": 1717245886944,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17545076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16847063.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16937641.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "d245790c5f3b26d675ddca05585d84364f4454f4",
          "message": "format\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T21:31:31+08:00",
          "tree_id": "795b6e4e084210be5326d56d92ae6c758cbdc5f4",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/d245790c5f3b26d675ddca05585d84364f4454f4"
        },
        "date": 1717248832892,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 17002887.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 16275364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 16414135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "6a756727a97441ffff86487db1165b6c7e02d787",
          "message": "Update README\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-06-01T21:36:42+08:00",
          "tree_id": "fae136515a245054fa88d11492de8578d5c9bd75",
          "url": "https://github.com/SciFracX/FractionalDiffEq.jl/commit/6a756727a97441ffff86487db1165b6c7e02d787"
        },
        "date": 1717249145010,
        "tool": "julia",
        "benches": [
          {
            "name": "FLMM/Trapezoid",
            "value": 18096207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22963808\nallocs=265935\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/NewtonGregory",
            "value": 17524745.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22419184\nallocs=264978\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "FLMM/BDF",
            "value": 17460771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22336672\nallocs=264966\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}