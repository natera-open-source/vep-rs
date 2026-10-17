(function pageToc() {
    'use strict';

    const content = document.getElementById('mdbook-content');
    const main = content ? content.querySelector('main') : null;
    if (!main) {
        return;
    }

    let headings = Array.from(main.querySelectorAll('h2[id], h3[id]'));
    if (headings.length < 3) {
        return;
    }
    if (headings.length > 30) {
        const sections = headings.filter(heading => heading.tagName === 'H2');
        if (sections.length >= 3) {
            headings = sections;
        }
    }

    const list = document.createElement('ol');
    const links = new Map();
    headings.forEach(heading => {
        const item = document.createElement('li');
        item.className = 'pagetoc-' + heading.tagName.toLowerCase();
        const link = document.createElement('a');
        link.href = '#' + heading.id;
        // Copies the heading anchor's nodes, so a flag in its <code> keeps the
        // no-wrap rule that plain text would lose. A link written inside a
        // heading closes the anchor while the page is parsed, leaving it
        // empty; that heading contributes its text only.
        const anchor = heading.querySelector(':scope > a.header');
        if (anchor && anchor.textContent.trim() && !anchor.querySelector('a')) {
            anchor.childNodes.forEach(node => link.appendChild(node.cloneNode(true)));
        } else {
            link.textContent = heading.textContent.trim();
        }
        item.appendChild(link);
        list.appendChild(item);
        links.set(heading, link);
    });

    const nav = document.createElement('nav');
    nav.className = 'pagetoc';
    nav.setAttribute('aria-label', 'On this page');
    const wide = window.matchMedia('(min-width: 1280px)');
    const pageTitle = main.querySelector(':scope > h1');

    function render() {
        nav.textContent = '';
        if (wide.matches) {
            const title = document.createElement('p');
            title.className = 'pagetoc-title';
            title.textContent = 'On this page';
            nav.append(title, list);
            content.insertBefore(nav, main);
        } else {
            const details = document.createElement('details');
            const summary = document.createElement('summary');
            summary.textContent = 'On this page';
            details.append(summary, list);
            nav.append(details);
            if (pageTitle) {
                pageTitle.after(nav);
            } else {
                main.prepend(nav);
            }
        }
    }

    render();
    wide.addEventListener('change', render);

    const menuBarHeight = parseInt(
        getComputedStyle(document.documentElement).getPropertyValue('--menu-bar-height'), 10) || 50;
    let active = null;

    function keepVisible(link) {
        if (!wide.matches) {
            return;
        }
        const linkBox = link.getBoundingClientRect();
        const navBox = nav.getBoundingClientRect();
        if (linkBox.top < navBox.top) {
            nav.scrollTop -= navBox.top - linkBox.top;
        } else if (linkBox.bottom > navBox.bottom) {
            nav.scrollTop += linkBox.bottom - navBox.bottom;
        }
    }

    function setActive(heading) {
        if (heading === active) {
            return;
        }
        if (active) {
            links.get(active).classList.remove('active');
            links.get(active).removeAttribute('aria-current');
        }
        active = heading;
        if (heading) {
            const link = links.get(heading);
            link.classList.add('active');
            link.setAttribute('aria-current', 'location');
            keepVisible(link);
        }
    }

    const visible = new Set();

    function pick() {
        let chosen = headings.find(heading => visible.has(heading));
        if (!chosen) {
            const limit = menuBarHeight + 8;
            for (const heading of headings) {
                if (heading.getBoundingClientRect().top < limit) {
                    chosen = heading;
                } else {
                    break;
                }
            }
        }
        // At the foot of the page the last sections may never reach the band
        // the observer watches.
        const root = document.documentElement;
        if (window.innerHeight + window.scrollY >= root.scrollHeight - 2) {
            const shown = headings.filter(heading => heading.getBoundingClientRect().top < window.innerHeight);
            if (shown.length) {
                chosen = shown[shown.length - 1];
            }
        }
        setActive(chosen || null);
    }

    // The observer reports only headings that cross its band between two
    // frames, so a jump past the band (Page Down, End, a restored position)
    // reaches the pick through the scroll listener.
    let queued = false;
    function schedule() {
        if (!queued) {
            queued = true;
            requestAnimationFrame(() => {
                queued = false;
                pick();
            });
        }
    }
    window.addEventListener('scroll', schedule, { passive: true });
    window.addEventListener('resize', schedule);

    if ('IntersectionObserver' in window) {
        const observer = new IntersectionObserver(entries => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    visible.add(entry.target);
                } else {
                    visible.delete(entry.target);
                }
            });
            pick();
        }, { rootMargin: '-' + menuBarHeight + 'px 0px -40% 0px', threshold: 0 });
        headings.forEach(heading => observer.observe(heading));
    }
    schedule();
})();
