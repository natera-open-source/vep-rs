/* global path_to_root */
(function vepRsTheme() {
    'use strict';

    const html = document.documentElement;
    const content = document.getElementById('mdbook-content');
    const main = content ? content.querySelector('main') : null;

    function stored(key) {
        try {
            return localStorage.getItem(key);
        } catch {
            return null;
        }
    }

    (function skipLink() {
        if (!content) {
            return;
        }
        content.setAttribute('tabindex', '-1');
        const link = document.createElement('a');
        link.className = 'skip-link';
        // The 404 page carries a <base href> of the site root, against which a
        // bare fragment would resolve to the home page.
        link.href = window.location.pathname + window.location.search + '#mdbook-content';
        link.textContent = 'Skip to content';
        document.body.insertBefore(link, document.body.firstChild);

        // toc.js calls scrollIntoView on the active sidebar entry, and Chromium
        // moves the sequential focus navigation starting point to that element,
        // so the first Tab would land after it. Focusing and blurring an inert
        // element placed before the skip link puts the starting point back.
        if (document.activeElement === document.body && !window.location.hash) {
            const anchor = document.createElement('div');
            anchor.tabIndex = -1;
            anchor.className = 'sr-only';
            anchor.setAttribute('aria-hidden', 'true');
            document.body.insertBefore(anchor, link);
            anchor.focus({ preventScroll: true });
            anchor.blur();
        }
    })();

    (function sidebar() {
        const checkbox = document.getElementById('mdbook-sidebar-toggle-anchor');
        const toggle = document.getElementById('mdbook-sidebar-toggle');
        const nav = document.getElementById('mdbook-sidebar');
        if (!checkbox || !toggle || !nav) {
            return;
        }
        // The same boundary as the overlay media query in vep-rs.css; a media
        // query counts the viewport with its scrollbar, as the stylesheet does.
        const pushLayout = window.matchMedia('(min-width: 1024px)');

        // Mirrors book.js's show and hide paths without writing the choice to
        // localStorage, so a default or a tap on the dimming layer is never
        // remembered as the reader's preference.
        function setOpen(open) {
            if (open) {
                nav.style.display = '';
                // Forces a reflow between the display change and the transform
                // transition, which Safari otherwise skips.
                nav.offsetHeight;
            }
            checkbox.checked = open;
            html.classList.toggle('sidebar-visible', open);
            nav.querySelectorAll('a').forEach(link => {
                link.setAttribute('tabIndex', open ? 0 : -1);
            });
            toggle.setAttribute('aria-expanded', String(open));
            nav.setAttribute('aria-hidden', String(!open));
        }

        if (pushLayout.matches && !checkbox.checked && stored('mdbook-sidebar') !== 'hidden') {
            // The stock inline script hid the sidebar before first paint; the
            // sidebar-resizing class holds the stock transitions off while the
            // page wrapper takes its margin, and a forced layout commits that
            // state before the class goes.
            const wrapper = document.getElementById('mdbook-page-wrapper');
            html.classList.add('sidebar-resizing');
            setOpen(true);
            if (wrapper) {
                wrapper.offsetWidth;
            }
            html.classList.remove('sidebar-resizing');
        }

        toggle.setAttribute('tabindex', '0');
        toggle.setAttribute('role', 'button');
        toggle.addEventListener('keydown', event => {
            if (event.key === 'Enter' || event.key === ' ') {
                event.preventDefault();
                toggle.click();
            }
        });

        const dim = document.createElement('div');
        dim.className = 'sidebar-dim';
        dim.setAttribute('aria-hidden', 'true');
        nav.parentNode.insertBefore(dim, nav.nextSibling);
        dim.addEventListener('click', () => setOpen(false));

        document.addEventListener('keydown', event => {
            const target = event.target;
            const typing = target && ['INPUT', 'TEXTAREA', 'SELECT'].includes(target.tagName);
            if (event.key === 'Escape' && checkbox.checked && !typing && !pushLayout.matches) {
                setOpen(false);
                toggle.focus();
            }
        });
    })();

    if (!main) {
        return;
    }

    (function longInlineCode() {
        const all = Array.from(content.querySelectorAll(':not(pre) > code'));
        // Code in a table cell is one unbreakable run, spaces included, with
        // a cut after a comma or semicolon that ends a piece of 16 characters
        // or more; scrollEdges activates those cuts for a table wider than
        // its wrapper, so a list of terms wraps and the table narrows, and a
        // span that fits stays on one line. Prose code, and code in the page
        // list of headings, wraps when it is wider than its block: between
        // its tokens and after a separator, and inside a token when the token
        // alone is wider than the block.
        const cellCode = all.filter(code => code.closest('.table-wrapper'));
        const spans = all.filter(code => !code.closest('.table-wrapper'));
        const listCuts = /(?<=[,;])/;
        const pathCuts = /(?<=[\/,;:=])/;
        const wideCuts = /(?<=[_.])/;

        function runOf(text) {
            const run = document.createElement('span');
            run.className = 'code-run';
            run.textContent = text;
            return run;
        }

        function placeRuns(code, pieces) {
            const runs = document.createDocumentFragment();
            pieces.forEach((piece, index) => {
                if (index && piece !== ' ') {
                    runs.append(document.createElement('wbr'));
                }
                runs.append(piece === ' ' ? ' ' : runOf(piece));
            });
            code.textContent = '';
            code.append(runs);
            code.setAttribute('data-runs', '');
        }

        // Splits a span into runs that never break inside: its
        // whitespace-separated tokens, each cut after a separator. A quoted
        // token keeps its slashes whole, so a sed or awk program breaks at a
        // space.
        function splitRuns(code) {
            if (code.hasAttribute('data-runs') || code.children.length) {
                return;
            }
            const pieces = [];
            code.textContent.split(/(\s+)/).forEach(part => {
                if (!part) {
                    return;
                }
                if (/^\s+$/.test(part)) {
                    pieces.push(' ');
                    return;
                }
                pieces.push(...part.split(/["']/.test(part) ? listCuts : pathCuts));
            });
            placeRuns(code, pieces);
        }

        function splitCellRuns(code) {
            if (code.hasAttribute('data-runs') || code.children.length) {
                return;
            }
            const pieces = [];
            let current = '';
            code.textContent.split(listCuts).forEach(piece => {
                current += piece;
                if (current.trim().length >= 16 && /[,;]$/.test(current)) {
                    pieces.push(current);
                    current = '';
                }
            });
            if (current) {
                pieces.push(current);
            }
            placeRuns(code, pieces);
        }

        cellCode.forEach(splitCellRuns);
        if (!spans.length) {
            return;
        }

        const blocks = new Map(spans.map(code => {
            let block = code.parentElement;
            while (block !== main && ['inline', 'contents'].includes(getComputedStyle(block).display)) {
                block = block.parentElement;
            }
            return [code, block];
        }));

        function blockWidth(code) {
            const block = blocks.get(code);
            const style = getComputedStyle(block);
            return block.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight);
        }

        // Only a span wider than its block wraps, between its runs; a run that
        // is itself wider than the block may break at any character.
        function fit() {
            spans.forEach(code => code.classList.remove('code-wrap'));
            main.querySelectorAll('.code-run-wide').forEach(run => run.classList.remove('code-run-wide'));
            const wide = spans.filter(code => code.getBoundingClientRect().width > blockWidth(code) + 0.5);
            wide.forEach(code => {
                splitRuns(code);
                code.classList.add('code-wrap');
            });
            wide.forEach(code => {
                // A run at the start or end of a line also carries the span's
                // side padding and border.
                const style = getComputedStyle(code);
                const sides = ['paddingLeft', 'paddingRight', 'borderLeftWidth', 'borderRightWidth']
                    .reduce((sum, key) => sum + (parseFloat(style[key]) || 0), 0);
                const limit = blockWidth(code) - sides;
                code.querySelectorAll('.code-run').forEach(run => {
                    if (run.getBoundingClientRect().width <= limit) {
                        return;
                    }
                    const pieces = run.textContent.split(wideCuts);
                    if (pieces.length === 1) {
                        run.classList.add('code-run-wide');
                        return;
                    }
                    const parts = document.createDocumentFragment();
                    pieces.forEach((piece, index) => {
                        if (index) {
                            parts.append(document.createElement('wbr'));
                        }
                        parts.append(runOf(piece));
                    });
                    run.replaceWith(parts);
                });
                code.querySelectorAll('.code-run').forEach(run => {
                    if (run.getBoundingClientRect().width > limit) {
                        run.classList.add('code-run-wide');
                    }
                });
            });
        }

        let queued = false;
        function refit() {
            if (!queued) {
                queued = true;
                requestAnimationFrame(() => {
                    queued = false;
                    fit();
                });
            }
        }

        fit();
        if ('ResizeObserver' in window) {
            new ResizeObserver(refit).observe(main);
        } else {
            window.addEventListener('resize', refit);
        }
        if (document.fonts && document.fonts.ready) {
            document.fonts.ready.then(refit);
        }
    })();

    (function scrollEdges() {
        function watch(scroller, host, onResize) {
            function update() {
                const max = scroller.scrollWidth - scroller.clientWidth;
                const states = [];
                if (max > 1) {
                    if (scroller.scrollLeft > 1) {
                        states.push('left');
                    }
                    if (scroller.scrollLeft < max - 1) {
                        states.push('right');
                    }
                }
                if (states.length) {
                    host.setAttribute('data-scroll', states.join(' '));
                    scroller.setAttribute('tabindex', '0');
                } else {
                    host.removeAttribute('data-scroll');
                    scroller.removeAttribute('tabindex');
                }
            }
            function resized() {
                if (onResize) {
                    onResize();
                }
                update();
            }
            scroller.addEventListener('scroll', update, { passive: true });
            if ('ResizeObserver' in window) {
                new ResizeObserver(resized).observe(scroller);
            } else {
                window.addEventListener('resize', resized);
            }
            resized();
            // Web fonts change a table's width without changing its wrapper's.
            if (document.fonts && document.fonts.ready) {
                document.fonts.ready.then(resized);
            }
        }

        // Below 480 px vep-rs.css gives a long cell a minimum width, so a table
        // too wide for the screen scrolls instead of squeezing its prose.
        const narrow = window.matchMedia('(max-width: 479px)');
        main.querySelectorAll('.table-wrapper td:not(:first-child)').forEach(cell => {
            const text = cell.textContent.trim();
            if (text.length > 40) {
                cell.classList.add('cell-long');
            } else if (text.length <= 20) {
                cell.classList.add(/\s/.test(text) ? 'cell-brief' : 'cell-short');
            }
        });
        // A column, parameter or property table stacks its rows on a phone;
        // each cell learns its column header for the stacked label.
        main.querySelectorAll('.table-wrapper table').forEach(table => {
            const heads = Array.from(table.querySelectorAll('thead th')).map(th => th.textContent.trim());
            if (heads.length < 2 || heads.length > 3
                || !['Column', 'Parameter', 'Property', 'Input', 'Plugin'].includes(heads[0])) {
                return;
            }
            table.classList.add('stack-rows');
            table.querySelectorAll('tbody tr').forEach(row => {
                Array.from(row.children).forEach((cell, index) => {
                    if (index && heads[index]) {
                        cell.setAttribute('data-label', heads[index]);
                    }
                });
            });
        });
        main.querySelectorAll('.table-wrapper thead').forEach(head => {
            if (!head.textContent.trim()) {
                head.classList.add('empty-head');
            }
        });
        main.querySelectorAll('.table-wrapper table').forEach(table => {
            const labels = Array.from(table.querySelectorAll('tbody tr > :first-child'));
            if (labels.length && labels.every(cell => cell.textContent.trim().length <= 24)) {
                table.classList.add('short-labels');
            }
        });

        main.querySelectorAll('.table-wrapper').forEach(wrapper => {
            const host = document.createElement('div');
            host.className = 'table-scroll';
            wrapper.parentNode.insertBefore(host, wrapper);
            host.appendChild(wrapper);
            // A blank header is hidden, so its cell has no width.
            const firstCell = wrapper.querySelector('thead:not(.empty-head) th:first-child, td:first-child');
            const longCell = wrapper.querySelector('td.cell-long');
            const values = Array.from(wrapper.querySelectorAll('tr > :not(:first-child)'));
            const table = wrapper.querySelector('table');
            const cutCode = Array.from(wrapper.querySelectorAll('code[data-runs]'))
                .filter(code => code.querySelector('wbr'));
            // A value that is one token with no cut point, such as a repeat-
            // expanded allele, may break at any character once it is wider
            // than the room beside the first column, or, for a token of 40
            // characters or more, wider than that room shared equally by the
            // value columns while the table still overflows.
            const soloCode = Array.from(wrapper.querySelectorAll('tr > :not(:first-child) code'))
                .filter(code => !/[\s,;]/.test(code.textContent) && !code.querySelector('wbr'));
            const firstRow = wrapper.querySelector('tr');
            const valueColumns = firstRow ? Math.max(1, firstRow.children.length - 1) : 1;

            function padding(cell) {
                const style = getComputedStyle(cell);
                return parseFloat(style.paddingLeft) + parseFloat(style.paddingRight);
            }

            function nearlyFits() {
                return wrapper.scrollWidth <= wrapper.clientWidth * 1.1;
            }

            function cellRoom(code) {
                const cell = code.closest('td, th');
                return cell.clientWidth - padding(cell);
            }

            function overflows() {
                return wrapper.scrollWidth > wrapper.clientWidth + 1;
            }

            watch(wrapper, host, () => {
                // Every span is measured with its classes removed, so its own
                // wrapping never feeds its test; the overflow checks read the
                // live layout, so a class already applied counts towards
                // fitting the table. Cell code breaks at its
                // cuts in an overflowing table, and there for a span wider
                // than the cell its column has; a span that fits moves to the
                // next line whole. In a stacked table each cell is its own
                // block to fit. A table that overflows with its headers wrapped
                // and is within a tenth of fitting gets tighter cell padding,
                // so a few pixels of scroll never squeeze the other columns.
                wrapper.classList.remove('wrap-heads', 'tight');
                cutCode.forEach(code => code.classList.remove('code-wrap'));
                soloCode.forEach(code => code.classList.remove('code-wide'));
                const widths = cutCode.map(code => code.getBoundingClientRect().width);
                const stacked = narrow.matches && table && table.classList.contains('stack-rows');
                const labelWidth = firstCell ? firstCell.offsetWidth : 0;
                // Widest first: once the widest token breaks and the table
                // fits, overflows() spares the narrower ones.
                const solo = soloCode.map(code => ({ code, width: code.getBoundingClientRect().width }))
                    .sort((a, b) => b.width - a.width);
                // Below 480 px an unstacked table scrolls whole, so nothing wraps.
                if (stacked || !narrow.matches) {
                    solo.forEach(entry => {
                        const cell = entry.code.closest('td, th');
                        const room = (stacked ? cell.clientWidth : wrapper.clientWidth - labelWidth) - padding(cell);
                        if (entry.width > room) {
                            entry.code.classList.add('code-wide');
                        }
                    });
                }
                if (!narrow.matches) {
                    solo.forEach(entry => {
                        const cell = entry.code.closest('td, th');
                        const share = (wrapper.clientWidth - labelWidth) / valueColumns - padding(cell);
                        if (entry.code.textContent.length >= 40 && entry.width > share && overflows()) {
                            entry.code.classList.add('code-wide');
                        }
                    });
                }
                if (stacked) {
                    cutCode.forEach((code, index) => {
                        if (widths[index] > cellRoom(code) + 0.5) {
                            code.classList.add('code-wrap');
                        }
                    });
                    // A stacked cell that still scrolls is reachable by keyboard.
                    table.querySelectorAll('td').forEach(cell => {
                        if (cell.scrollWidth > cell.clientWidth + 1) {
                            cell.setAttribute('tabindex', '0');
                        } else {
                            cell.removeAttribute('tabindex');
                        }
                    });
                } else if (cutCode.length && overflows()) {
                    cutCode.forEach(code => code.classList.add('code-wrap'));
                    cutCode.forEach((code, index) => {
                        if (widths[index] <= cellRoom(code) + 0.5) {
                            code.classList.remove('code-wrap');
                        }
                    });
                }
                if (!narrow.matches) {
                    if (overflows()) {
                        wrapper.classList.add('wrap-heads');
                    }
                    if (overflows() && nearlyFits()) {
                        wrapper.classList.add('tight');
                    }
                    if (table) {
                        table.querySelectorAll('td[tabindex]').forEach(cell => cell.removeAttribute('tabindex'));
                    }
                }
                if (!firstCell) {
                    return;
                }
                const first = firstCell.offsetWidth;
                const room = wrapper.clientWidth - first;
                let sticky = false;
                if (narrow.matches) {
                    // No pinned column on a phone; a long cell keeps a floor of
                    // 12em so the table scrolls instead of squeezing its prose.
                    if (longCell) {
                        const style = getComputedStyle(longCell);
                        const floor = 12 * parseFloat(style.fontSize)
                            + parseFloat(style.paddingLeft) + parseFloat(style.paddingRight);
                        wrapper.style.setProperty('--cell-floor', Math.floor(floor) + 'px');
                    }
                } else {
                    // The first column is pinned when it leaves most of the
                    // wrapper free and every other cell fits beside it, so no
                    // value is hidden at every scroll position.
                    const widest = values.reduce((width, cell) => Math.max(width, cell.getBoundingClientRect().width), 0);
                    sticky = first <= 0.6 * wrapper.clientWidth && widest <= room;
                    wrapper.style.removeProperty('--cell-floor');
                }
                wrapper.classList.toggle('no-sticky', !sticky);
                host.style.setProperty('--sticky-width', (sticky ? first : 0) + 'px');
            });
        });

        main.querySelectorAll('pre').forEach(pre => {
            const code = pre.querySelector(':scope > code');
            if (code) {
                watch(code, pre);
            }
        });
    })();

    (function menuTitle() {
        // The book title in the menu bar is an h1 in the template; the page's
        // own title is the one heading of that level.
        const title = document.querySelector('h1.menu-title');
        if (!title) {
            return;
        }
        const label = document.createElement('p');
        label.className = title.className;
        label.textContent = title.textContent;
        label.addEventListener('click', () => {
            document.scrollingElement.scrollTo({ top: 0, behavior: 'smooth' });
        });
        title.replaceWith(label);
    })();

    (function editLink() {
        // The edit template only knows the chapter file, a stub whose text the
        // preprocessor takes from the repository file.
        const edit = document.querySelector('#mdbook-menu-bar a[rel="edit"]');
        if (!edit) {
            return;
        }
        const sources = { 'docs/changelog.md': 'CHANGELOG.md', 'docs/contributing.md': 'CONTRIBUTING.md' };
        Object.keys(sources).forEach(stub => {
            if (edit.href.endsWith('/' + stub)) {
                edit.href = edit.href.slice(0, -stub.length) + sources[stub];
            }
        });
    })();

    (function channels() {
        const root = typeof path_to_root === 'string' ? path_to_root : '';
        const parentSiteHref = root + '../';
        const channelMeta = document.querySelector('meta[name="vep-rs-channel"]');
        const channel = channelMeta ? channelMeta.getAttribute('content') : null;

        function current(channels) {
            const here = window.location.pathname;
            let best = null;
            channels.forEach(entry => {
                if (here.startsWith(entry.path) && (!best || entry.path.length > best.path.length)) {
                    best = entry;
                }
            });
            return best;
        }

        function banner(releaseHref) {
            if (channel !== 'main' || !content || content.querySelector('.vep-banner')) {
                return;
            }
            const note = document.createElement('div');
            note.className = 'vep-banner';
            note.setAttribute('role', 'note');
            const link = document.createElement('a');
            link.href = releaseHref;
            link.textContent = 'the release documentation';
            note.append('Development version from the main branch; see ', link, '.');
            content.insertBefore(note, content.firstChild);
        }

        function selector(channels, release) {
            const byPath = current(channels);
            // A site served at a path versions.json does not list still knows
            // its channel: the preview's meta names it, a release page has none.
            const here = byPath || (channel
                ? channels.find(entry => entry.name.split(' ')[0] === channel)
                : release) || null;
            const buttons = document.querySelector('#mdbook-menu-bar .right-buttons');
            if (channels.length < 2 || !here || !buttons) {
                return;
            }
            const label = document.createElement('label');
            label.className = 'vep-version';
            const text = document.createElement('span');
            text.className = 'sr-only';
            text.textContent = 'Documentation version';
            const select = document.createElement('select');
            channels.forEach(entry => {
                const option = document.createElement('option');
                option.value = entry.path;
                option.textContent = entry.name;
                option.selected = entry === here;
                select.appendChild(option);
            });
            select.addEventListener('change', () => {
                if (!byPath) {
                    window.location.href = select.value;
                    return;
                }
                const page = window.location.pathname.slice(here.path.length);
                const target = select.value + page + window.location.hash;
                fetch(select.value + page, { method: 'HEAD' })
                    .then(response => {
                        window.location.href = response.ok ? target : select.value;
                    })
                    .catch(() => {
                        window.location.href = select.value;
                    });
            });
            label.append(text, select);
            buttons.insertBefore(label, buttons.firstChild);
            document.getElementById('mdbook-menu-bar').classList.add('has-version');
        }

        fetch(root + 'versions.json')
            .then(response => (response.ok ? response.json() : null))
            .then(data => {
                if (!data) {
                    banner(parentSiteHref);
                    return;
                }
                const channels = Array.isArray(data.channels) ? data.channels : [];
                const release = channels.find(entry => entry.name === data.release);
                banner(release ? release.path : parentSiteHref);
                selector(channels, release);
            })
            .catch(() => banner(parentSiteHref));
    })();
})();
