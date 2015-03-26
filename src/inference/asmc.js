////////////////////////////////////////////////////////////////////
// Asynchronous Anytime SMC.
// http://arxiv.org/abs/1407.2864
// bufferSize: how many particles to keep in queue
// numParticles: total number of particles to run

'use strict';

var _ = require('underscore');
var util = require('../util.js');
var erp = require('../erp.js');

module.exports = function(env) {

  function copyOneParticle(particle) {
    return {
      continuation: particle.continuation,
      weight: particle.weight,
      finalWeight: particle.finalWeight,
      completed: particle.completed,
      factorIndex: particle.factorIndex,
      value: particle.value,
      numChildrenToSpawn: 1,
      multiplicity: particle.multiplicity,
      store: _.clone(particle.store)
    };
  }

  function initParticle(s, cont) {
    return {
      continuation: cont,
      weight: 0,
      finalWeight: 0,
      completed: false,
      factorIndex: undefined,
      value: undefined,
      numChildrenToSpawn: 0,
      multiplicity: 1,
      store: _.clone(s)
    };
  }

  function aSMC(s, k, a, wpplFn, numParticles, bufferSize) {

    this.buffer = [];
    this.bufferSize = bufferSize; // \rho
    this.numParticles = 0;

    this.obsWeights = {};
    this.exitedParticles = [];

    this.exitK = function(s) {return wpplFn(s, env.exit, a);};
    this.store = s;
    console.log("this.store = ", this.store);

    // Create initial particles
    var initNumP = Math.floor(bufferSize * 3 / 5) // \rho_0 = 3/5 * \rho
    for (var i = 0; i < initNumP; i++) {
      this.buffer.push(initParticle(s, this.exitK));
    }

    // Move old coroutine out of the way and install this as the current
    // handler.
    this.k = k;
    this.oldCoroutine = env.coroutine;
    env.coroutine = this;
    this.oldStore = _.clone(s); // will be reinstated at the end

    return this.control(numParticles)
  };

  // aSMC.prototype.run = function(numParticles) {
  //   // Run first particle
  //   return this.control(numParticles)
  // };

  aSMC.prototype.control = function(numP) {
    // allows for continuation of shooting particles
    this.numParticles = (numP == undefined) ? this.numParticles : this.numParticles + numP;

    // make a choice of whether you want to launch a new particle or continue
    // an existing one
    var i = Math.floor((this.buffer.length + 1) * Math.random())
    var launchP = i == this.buffer.length ? initParticle(this.store, this.exitK) : this.buffer[i];
    var p;

    if (launchP.numChildrenToSpawn > 1) { // launch one and reduce numChildren
      p = copyOneParticle(launchP);
      launchP.numChildrenToSpawn -= 1;
    } else {           // if initial particle or particle with 1 child to spawn
      p = launchP;
      this.buffer = _.without(this.buffer, p) // deque
    }
    this.activeParticle = p;
    return p.continuation(p.store)
  };

  aSMC.prototype.sample = function(s, cc, a, erp, params) {
    return cc(s, erp.sample(params));
  };

  aSMC.prototype.factor = function(s, cc, a, score) {
    this.activeParticle.weight += score;
    this.activeParticle.continuation = cc;
    this.activeParticle.store = s;
    var fi = this.activeParticle.factorIndex
    var newFI = fi == undefined ? 0 : fi + 1;
    this.activeParticle.factorIndex = newFI;

    var lk = this.obsWeights[newFI];
    if (lk == undefined) {      // 1st particle at observation
      var det = {
        w: this.activeParticle.weight,
        wbar: this.activeParticle.weight,
        mnk:  1,
        vnk: this.activeParticle.weight
      };
      this.obsWeights[newFI] = [det];
      this.activeParticle.finalWeight = this.activeParticle.weight;
    } else {                    // 2nd or more particle at observation
      var currMultiplicity = this.activeParticle.multiplicity;
      var currWeight = this.activeParticle.weight;
      var ldenom = Math.log(lk.length + currMultiplicity); // k - 1 + Ckn : log
      var prevWBar = lk[lk.length-1].wbar;
      var wbar = util.logsumexp([Math.log(lk.length) - ldenom + prevWBar,
                                 Math.log(currMultiplicity - ldenom + currWeight)])
      var logRatio = currWeight - wbar; // Rnk : log
      var numChildrenAndWeight = [];
      if (logRatio < 0) {
        numChildrenAndWeight = Math.log(Math.random()) < logRatio ?
          [1, wbar] :
          [0, -Infinity];
      } else {
        var totalChildren = 0;
        for (var v = 0; v < lk.length; v++) totalChildren += lk[v].mnk // \sum M^k_n
        var minK = Math.min(this.bufferSize, lk.length);               // min(K_0, k-1)
        var rnk = Math.exp(logRatio);
        var clampedRnk = totalChildren <= minK ? Math.ceil(rnk) : Math.floor(rnk);
        numChildrenAndWeight = [clampedRnk, currWeight - Math.log(clampedRnk)];
      }
      var det2 = {
        w: currWeight,
        wbar: wbar,
        mnk:  numChildrenAndWeight[0],
        vnk: numChildrenAndWeight[1]
      };
      this.obsWeights[newFI] = lk.concat([det2])
      // console.log("this.obsWeights = ", this.obsWeights);
      if (numChildrenAndWeight[0] > 0) {            // there are children
        if (this.buffer.length < this.bufferSize) { // buffer can be added to
          this.activeParticle.numChildrenToSpawn = numChildrenAndWeight[0];
          this.activeParticle.weight = numChildrenAndWeight[1];
        } else {                  // buffer full, update multiplicty
          this.activeParticle.multiplicity *= numChildrenAndWeight[0];
          this.activeParticle.numChildrenToSpawn = 1;
          this.activeParticle.weight = numChildrenAndWeight[1];
        }
        // when at last factor, this weight will be the correct weight for the particle
        this.activeParticle.finalWeight = Math.log(this.activeParticle.currMultiplicity)
          + this.activeParticle.weight
          + score
        this.buffer.push(this.activeParticle); // push into buffer
      }
      // if there are no children, active particle automatically dropped
    }
    return this.control()              // return to control
  };

  aSMC.prototype.exit = function(s, retval) {
    // needswork: if numParticles number of particles have exited, then package
    // up the continuation and return to wppl
    this.activeParticle.value = retval;
    console.log("this.activeParticle.value = ", this.activeParticle.value);
    this.activeParticle.completed = true;

    this.activeParticle.weight = this.activeParticle.finalWeight;
    this.exitedParticles.push(this.activeParticle)
    console.log("this.exitedParticles = ", this.exitedParticles);

    if (this.exitedParticles.length < this.numParticles) {
      return this.control()
    } else {
      // Compute marginal distribution from (unweighted) particles
      var hist = {};
      _.each(
        this.exitedParticles,
        function(particle) {
          var k = JSON.stringify(particle.value);
          if (hist[k] === undefined) {
            hist[k] = { prob: 0, val: particle.value };
          }
          hist[k].prob += 1;
        });
        var dist = erp.makeMarginalERP(hist);

      var nc = util.logsumexp(_.map(this.exitedParticles, function(p) {return p.weight}));
      // Save estimated normalization constant in erp (average particle weight)
      dist.normalizationConstant = nc - Math.log(this.numParticles);
      dist.continue = function(numP) {this.control(numP)};
      // Reinstate previous coroutine:
      env.coroutine = this.oldCoroutine;
      // Return from particle filter by calling original continuation:
      return this.k(this.oldStore, dist);
    }
  };

  function asyncPF(s, cc, a, wpplFn, numParticles, bufferSize) {
    return new aSMC(s, cc, a, wpplFn, numParticles, bufferSize);
  }

  return {
    aSMC: asyncPF
  };

};
