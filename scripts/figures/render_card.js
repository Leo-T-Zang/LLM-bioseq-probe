// Renders scripts/figures/social_card.html to figures/social_card.png (1200x630 at 2x).
// Needs playwright-core and a Chromium: CHROME=/path/to/chrome node scripts/figures/render_card.js
const path = require('path');
const { chromium } = require(process.env.PW || 'playwright-core');
(async () => {
  const b = await chromium.launch({ executablePath: process.env.CHROME, args: ['--no-sandbox'] });
  const p = await b.newPage({ viewport: { width: 1200, height: 630 }, deviceScaleFactor: 2 });
  await p.goto('file://' + path.resolve(__dirname, 'social_card.html'), { waitUntil: 'networkidle' });
  await p.evaluate(() => document.fonts.ready);
  await p.screenshot({ path: path.resolve(__dirname, '../../figures/social_card.png') });
  await b.close();
})();
