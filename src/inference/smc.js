////////////////////////////////////////////////////////////////////
// Particle filter with lightweight MH rejuvenation.
//
// Sequential importance re-sampling, which treats 'factor' calls as
// the synchronization / intermediate distribution points.
//
// After each factor particles are rejuvenated via lightweight MH.
// If numParticles==1 this amounts to MH with an (expensive) annealed init (returning only one sample)
// If rejuvSteps==0 this is a plain PF without any MH.

'use strict';

var _ = require('underscore');
var assert = require('assert');
var util = require('../util.js');
var erp = require('../erp.js');

module.exports = function(env) {

  var mh = require('./mh.js')(env);

  function copyPFRParticle(particle) {
    return {
      continuation: particle.continuation,
      weight: particle.weight,
      score: particle.score,
      trace: mh.deepCopyTrace(particle.trace),
      restrictedRegenFrom: particle.restrictedRegenFrom,
      completed: particle.completed,
      value: particle.value,
      store: _.clone(particle.store)
    };
  }

  function ParticleFilterRejuv(s, k, a, wpplFn, numParticles, rejuvSteps, restrict) {
    this.particles = [];
    this.particleIndex = 0;  // marks the active particle
    this.rejuvSteps = rejuvSteps;
    this.baseAddress = a;
    this.wpplFn = wpplFn;
    this.isParticleFilterRejuvCoroutine = true;
    this.restrict = restrict || false;

    // Create initial particles
    var exitK = function(s) {return wpplFn(s, env.exit, a);};
    for (var i = 0; i < numParticles; i++) {
      var particle = {
        continuation: exitK,
        weight: 0,
        score: 0,
        trace: [],
        restrictedRegenFrom: 0,
        completed: false,
        value: undefined,
        store: _.clone(s)
      };
      this.particles.push(particle);
    }

    // Move old coroutine out of the way and install this as current handler.
    this.k = k;
    this.oldCoroutine = env.coroutine;
    env.coroutine = this;
    this.oldStore = s; // will be reinstated at the end
  }

  ParticleFilterRejuv.prototype.run = function() {
    return this.activeParticle().continuation(this.activeParticle().store);
  };

  ParticleFilterRejuv.prototype.sample = function(s, cc, a, erp, params) {
    var val = erp.sample(params);
    var currScore = this.activeParticle().score;
    var choiceScore = erp.score(params, val);
    this.activeParticle().trace.push(
        mh.makeTraceEntry(_.clone(s), cc, a, erp, params, currScore, choiceScore, val, false)
    );
    this.activeParticle().score += choiceScore;
    return cc(s, val);
  };

  ParticleFilterRejuv.prototype.factor = function(s, cc, a, score) {
    // Update particle weight and score
    this.activeParticle().weight += score;
    this.activeParticle().score += score;
    this.activeParticle().continuation = cc;
    this.activeParticle().store = _.clone(s);

    if (this.allParticlesAdvanced()) {
      // Resample in proportion to weights
      this.resampleParticles();
      //rejuvenate each particle via MH
      return util.cpsForEach(
          function(particle, i, particles, nextK) {
            // make sure mhp coroutine doesn't escape:
            assert(env.coroutine.isParticleFilterRejuvCoroutine);
            return new MHP(function(p) {
              if (env.coroutine.restrict) p.restrictedRegenFrom = p.trace.length;
              particles[i] = p;
              return nextK();
            }, particle, this.baseAddress, a, this.wpplFn, this.rejuvSteps).run();
          }.bind(this),
          function() {
            // variable #factors: resampling can kill all continuing particles
            var fp = this.firstRunningParticleIndex();
            this.particleIndex = fp == -1 ? this.particles.length - 1 : fp;
            return this.activeParticle().continuation(this.activeParticle().store);
          }.bind(this),
          this.particles
      );
    } else {
      // Advance to the next particle
      this.particleIndex = this.nextRunningParticleIndex();
      return this.activeParticle().continuation(this.activeParticle().store);
    }
  };

  ParticleFilterRejuv.prototype.activeParticle = function() {
    return this.particles[this.particleIndex];
  };

  ParticleFilterRejuv.prototype.firstRunningParticleIndex = function() {
    return util.indexOfPred(this.particles, function(p) {return !p.completed});
  };

  ParticleFilterRejuv.prototype.nextRunningParticleIndex = function() {
    var ni = this.particleIndex + 1;
    var nxt = util.indexOfPred(this.particles, function(p) {return !p.completed}, ni);
    return (nxt >= 0 ?
            nxt :
            util.indexOfPred(this.particles, function(p) {return !p.completed}));
  };

  ParticleFilterRejuv.prototype.lastRunningParticleIndex = function() {
    return util.lastIndexOfPred(this.particles, function(p) {return !p.completed});
  };

  ParticleFilterRejuv.prototype.allParticlesAdvanced = function() {
    return this.particleIndex === this.lastRunningParticleIndex();
  };

  ParticleFilterRejuv.prototype.resampleParticles = function() {
    // Residual resampling following Liu 2008; p. 72, section 3.4.4
    var m = this.particles.length;
    var W = util.logsumexp(_.map(this.particles, function(p) {return p.weight;}));
    var avgW = W - Math.log(m);

    // Allow -Infinity case (for mh initialization, in particular with few particles)
    if (avgW == -Infinity) {
      console.warn('ParticleFilterRejuv: resampleParticles: all ' + m + ' particles have weight -Inf');
      return;
    }
    // Compute list of retained particles
    var retainedParticles = [];
    var newExpWeights = [];
    _.each(
        this.particles,
        function(particle) {
          var w = Math.exp(particle.weight - avgW);
          var nRetained = Math.floor(w);
          newExpWeights.push(w - nRetained);
          for (var i = 0; i < nRetained; i++) {
            retainedParticles.push(copyPFRParticle(particle));
          }});
    // Compute new particles
    var numNewParticles = m - retainedParticles.length;
    var j, newParticles = [];
    for (var i = 0; i < numNewParticles; i++) {
      j = erp.multinomialSample(newExpWeights);
      newParticles.push(copyPFRParticle(this.particles[j]));
    }
    // Particles after update: Retained + new particles
    this.particles = newParticles.concat(retainedParticles);
    // Reset all weights
    _.each(this.particles, function(particle) {particle.weight = avgW;});
  };

  ParticleFilterRejuv.prototype.exit = function(s, retval) {
    this.activeParticle().value = retval;
    this.activeParticle().completed = true;
    // this should be negative if there are no valid next particles
    var nextRunningParticleIndex = this.nextRunningParticleIndex();
    var allParticlesFinished = nextRunningParticleIndex < 0;

    // Wait for all particles to reach exit
    if (!allParticlesFinished) {
      this.particleIndex = nextRunningParticleIndex;
      return this.activeParticle().continuation(this.activeParticle().store);
    }

    // Final rejuvenation:
    var oldStore = this.oldStore;
    var hist = util.initHashMap();
    return util.cpsForEach(
        function(particle, i, particles, nextK) {
          // make sure mhp coroutine doesn't escape:
          assert(env.coroutine.isParticleFilterRejuvCoroutine);
          particle.restrictedRegenFrom = 0; // no restrictions on final run
          return new MHP(function(p) {particles[i] = p; return nextK();},
                       particle, this.baseAddress, undefined,
                       this.wpplFn, this.rejuvSteps, hist).run();
        }.bind(this),
        function() {
          var dist = erp.makeMarginalERP(hist);
          // Save estimated normalization constant in erp (average particle weight)
          dist.normalizationConstant = this.particles[0].weight;
          // Reinstate previous coroutine:
          var k = this.k;
          env.coroutine = this.oldCoroutine;
          // Return from particle filter by calling original continuation:
          return k(oldStore, dist);
        }.bind(this),
        this.particles
    );
  };

  ////// Lightweight MH on a particle

  function MHP(backToPF, particle, baseAddress, limitAddress, wpplFn, numIterations, hist) {
    this.oldStore = particle.store; // previous store at limitAddress
    this.val = particle.value;
    this.backToPF = backToPF;
    this.iterations = numIterations;
    this.totalIterations = numIterations;
    this.limitAddress = limitAddress;
    this.originalParticle = particle;
    this.hist = hist;
    this.regenFrom = 0;
    this.prefixTraceLength = particle.restrictedRegenFrom;
    this.remainingTrace = particle.trace.slice(this.prefixTraceLength);

    // Move PF coroutine out of the way and install this as current handler.
    this.oldCoroutine = env.coroutine;
    env.coroutine = this;
  }

  MHP.prototype.run = function() {
    this.sites = {};
    // if restricted, from first erp after last factor, else from beginning
    this.trace = this.remainingTrace;
    this.currScore = this.originalParticle.score;
    this.oldScore = -Infinity;
    this.fwdLP = 0;
    this.bwdLP = 0;

    // exit without doing rejuvenation when in restricted mode if:
    // a. the score coming in here is -Infinity or
    // b. there are no erps to sample between last and current factor
    if (this.oldCoroutine.restrict && (this.currScore == -Infinity || this.trace.length == 0)) {
      env.coroutine = this.oldCoroutine;
      return this.backToPF(this.originalParticle);
    };
    if (this.iterations === 0) {
      env.coroutine = this.oldCoroutine;
      return this.backToPF(this.originalParticle);
    } else { // FIXME: on final exit, will this end up calling the MH exit correctly?
      return this.propose(this.val, mh.deepCopyTrace);
    }
  };

  MHP.prototype.factor = function(s, k, a, sc) {
    this.currScore += sc;
    // exit if we've reached the fathest point of this particle
    return a === this.limitAddress ? env.exit(s) : k(s);
  };

  MHP.prototype.sample = function() {
    return mh.sampler.apply(this, arguments);
  };

  MHP.prototype.propose = function() {
    return mh.proposer.apply(this, arguments);
  };

  MHP.prototype.exit = function(s, val) {
    this.val = val;

    if (this.iterations > 0) {
      if (this.iterations === this.totalIterations && this.currScore === -Infinity)
        return this.run();      // when uninitialized and failing, do rejection-sampling
      this.iterations -= 1;

      // housekeeping - remove dead refs and recompute backward logprob
      if (this.oldTrace !== undefined && this.currScore !== -Infinity) {
        var i, reached = {};
        for (i = this.regenFrom; i < this.trace.length; i++)
          reached[this.trace[i].name] = this.trace[i];
        for (i = this.regenFrom; i < this.oldTrace.length; i++) {
          var v = this.oldTrace[i];
          if (reached[v.name] === undefined) {
            delete this.sites[v.name];
            this.bwdLP += v.choiceScore;
          }
        }
      }
      // Did we like this proposal?
      var acceptance = mh.acceptProb(this.currScore,
                                     this.oldScore,
                                     this.trace.length + this.prefixTraceLength,
                                     this.oldTrace === undefined ?
                                     0 :
                                     this.oldTrace.length + this.prefixTraceLength,
                                     this.bwdLP,
                                     this.fwdLP);
      if (Math.random() >= acceptance) { // If rejected, roll back trace, etc:
        this.sites = this.oldSites;
        this.trace = this.oldTrace;
        this.currScore = this.oldScore;
        this.val = this.oldVal;
      } else {
        this.oldStore = s;
      }
      // If final rejuvenation run, build hist from *all* MCMC steps
      if (this.hist !== undefined) {// Compute marginal distribution from (unweighted) particles
        var lk = this.hist.get(this.val);
        if (!lk) this.hist.set(this.val, 0);
        this.hist.set(this.val, this.hist.get(this.val) + 1);
      }

      return this.propose(this.val, mh.deepCopyTrace);
    } else {
      var newParticle = {
        continuation: this.originalParticle.continuation,
        weight: this.originalParticle.weight,
        score: this.currScore,
        // splice trace from first erp after last factor into original trace
        trace: this.originalParticle.trace.slice(0, this.prefixTraceLength).concat(this.trace),
        restrictedRegenFrom: this.originalParticle.restrictedRegenFrom,
        completed: this.originalParticle.completed,
        value: this.val,
        store: _.clone(this.oldStore) // use store from latest accepted proposal
      };
      // Reinstate previous coroutine and return by calling original continuation:
      env.coroutine = this.oldCoroutine;
      return this.backToPF(newParticle);
    }
  };

  function pfr(s, cc, a, wpplFn, numParticles, rejuvSteps, restrict) {
    return new ParticleFilterRejuv(s, cc, a, wpplFn, numParticles, rejuvSteps, restrict).run();
  }

  return {
    ParticleFilterRejuv: pfr
  };

};
