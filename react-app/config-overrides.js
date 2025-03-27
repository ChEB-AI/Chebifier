module.exports = function override(config, env) {
  console.log("React app rewired works!")
  config.resolve.fallback = {
    fs: false,
    path: false,
    polyfill: false,
  };
  return config;
};