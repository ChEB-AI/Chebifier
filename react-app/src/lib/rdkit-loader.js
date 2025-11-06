// src/lib/rdkit-loader.js
let loadPromise;

export function loadRDKit() {
  if (window.RDKit) return Promise.resolve(window.RDKit);
  if (loadPromise) return loadPromise;

  loadPromise = new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = process.env.PUBLIC_URL + '/RDKit_minimal.js';
    script.async = true;
    script.onload = () => {
      // initRDKitModule is defined by RDKit_minimal.js
      window.initRDKitModule()
        .then((RDKit) => {
          window.RDKit = RDKit; // cache globally
          resolve(RDKit);
        })
        .catch(reject);
    };
    script.onerror = reject;
    document.body.appendChild(script);
  });

  return loadPromise;
}